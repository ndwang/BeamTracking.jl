function _track!(
  i,
  coords::Coords,
  bunch::Bunch,
  ele::LineElement, 
  tm::SpaceChargeIntegration;
  kwargs...
)
  # Unpack the line element (type unstable)
  L = ele.L # Automatically calls deval (element-level get)
  ap = deval(ele.AlignmentParams)
  bp = deval(ele.BendParams)
  bm = deval(ele.BMultipoleParams)
  pp = deval(ele.PatchParams)
  dp = deval(ele.ApertureParams)
  sc = deval(ele.SpaceChargeParams)

  sc_calc_step = L / tm.num_sc_steps
  kc = fetch_kernels(i, coords, tm, bunch, sc_calc_step, ap, bp, bm, pp, dp, sc; kwargs...)
  init_efield_scratch(scp)

  for i_step in 1:tm.num_sc_steps
    SpaceChargeIntegrationTracking.sc_calc(sc, bunch)
    runkernels!(i, coords, kc; kwargs...)
  end
end

@inline function determine_element_type(
  kc,
  i, 
  coords,
  tm,
  bunch,
  L, 
  alignmentparams, 
  bendparams,
  bmultipoleparams,
  patchparams,
  apertureparams,
  spacechargeparams;
  kwargs...
)
  
  # Zero-length elements
  if L ≈ 0
    if isactive(bmultipoleparams)
      return @inline(thin_pure_bdipole(SplitIntegration(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bmultipoleparams))
    end
  end

  # Bend
  if isactive(bendparams)
    if !isactive(bmultipoleparams)
      # Bend no field or pure dipole
      return @inline(thick_bend_no_field(tm, bunch, bendparams, spacechargeparams, L))
    elseif 1 in bmultipoleparams.order && get_n_multipoles(bmultipoleparams) == 1
      # Pure dipole
      return @inline(thick_bend_pure_bdipole(tm, bunch, bendparams, first(bmultipoleparams), spacechargeparams, L))
    else
      # Bend with higher-order multipole components
      return @inline(bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L))
    end
  end

  # Drift
  if !isactive(bmultipoleparams)
    return @inline(drift(tm, bunch, spacechargeparams, L))
  end

  # Multipoles
  if 0 in bmultipoleparams.order
    # Solenoid
    if get_n_multipoles(bmultipoleparams) == 1
      return @inline(thick_pure_bsolenoid(tm, bunch, first(bmultipoleparams), spacechargeparams, L))
    else
      return @inline(thick_bsolenoid(tm, bunch, bmultipoleparams, spacechargeparams, L))
    end
  elseif 2 in bmultipoleparams.order 
    # Quadrupole
    return @inline(thick_bquadrupole(tm, bunch, bmultipoleparams, spacechargeparams, L))
  else
    # Steering or other higher-order multipole elements
    return @inline(thick_bmultipole(tm, bunch, bmultipoleparams, spacechargeparams, L))
  end
end

function fetch_kernels(
  i, 
  coords,
  tm,
  bunch,
  L, 
  alignmentparams, 
  bendparams,
  bmultipoleparams,
  patchparams,
  apertureparams,
  spacechargeparams;
  kwargs...
)
  kc = KernelChain(Val{4}())

  if isactive(patchparams)
    error("Tracking through a LineElement containing both PatchParams and SpaceChargeParams not currently defined")
  end

  if isactive(alignmentparams)
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, true)))
  end

  kc = push(kc, @inline(SpaceChargeIntegrationTracking.interpolate_field(i, coords, spacechargeparams)))
  kc = push(kc, @inline(determine_element_type(kc, i, coords, tm, bunch, L, alignmentparams, bendparams, bmultipoleparams, patchparams, apertureparams, spacechargeparams; kwargs...)))

  if isactive(alignmentparams)
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, false)))
  end

  return kc
end