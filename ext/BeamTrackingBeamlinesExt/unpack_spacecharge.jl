function _track!(
  i,
  coords::Coords,
  bunch::Bunch,
  ele::LineElement, 
  tm::SpaceCharge;
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

  sc_calc_step = L / sc.steps
  kc = fetch_kernels(i, coords, tm, bunch, sc_calc_step / 2, ap, bp, bm, pp, dp, sc; kwargs...)

  runkernels!(i, coords, kc; kwargs...)
  sc_calc() # Needs to be implemented
  runkernels!(i, coords, kc; kwargs...)

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
  kc = KernelChain(Val{1}())

  if isactive(patchparams)
    if isactive(bendparams) || isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both PatchParams and BendParams/BMultipoleParams not currently defined")
    else
      kc = push(kc, @inline(pure_patch(Exact(), bunch, patchparams, L)))
    end
  end

  # Zero-length elements
  if L ≈ 0
    if isactive(bmultipoleparams)
      kc = push(kc, @inline(thin_pure_bdipole(SplitIntegration(order=tm.order, num_steps=tm.num_steps, ds_step=tm.ds_step), bunch, bmultipoleparams)))
    end
  end

  # Bend
  if isactive(bendparams)
    if !isactive(bmultipoleparams)
      # Bend no field or pure dipole
      kc = push(kc, @inline(thick_bend_no_field(tm, bunch, bendparams, spacechargeparams, L)))
    elseif 1 in bmultipoleparams.order && get_n_multipoles(bmultipoleparams) == 1
      # Pure dipole
      kc = push(kc, @inline(thick_bend_pure_bdipole(tm, bunch, bendparams, first(bmultipoleparams), spacechargeparams, L)))
    else
      # Bend with higher-order multipole components
      kc = push(kc, @inline(bend_bmultipole(tm, bunch, bendparams, bmultipoleparams, L)))
    end
  end

  # Drift
  if !isactive(bmultipoleparams)
    kc = push(kc, @inline(drift(tm, bunch, spacechargeparams, L)))
  end

  # Multipoles
  if 0 in bmultipoleparams.order
    # Solenoid
    if get_n_multipoles(bmultipoleparams) == 1
      kc = push(kc, @inline(thick_pure_bsolenoid(tm, bunch, first(bmultipoleparams), spacechargeparams, L)))
    else
      kc = push(kc, @inline(thick_bsolenoid(tm, bunch, bmultipoleparams, spacechargeparams, L)))
    end
  elseif 2 in bmultipoleparams.order 
    # Quadrupole
    kc = push(kc, @inline(thick_bquadrupole(tm, bunch, bmultipoleparams, spacechargeparams, L)))
  else
    # Steering or other higher-order multipole elements
    kc = push(kc, @inline(thick_bmultipole(tm, bunch, bmultipoleparams, spacechargeparams, L)))
  end

  return kc
end