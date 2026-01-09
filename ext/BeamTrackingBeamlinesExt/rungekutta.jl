# =========== BYPASS UNPACKING FOR RUNGEKUTTA ============= #

"""
Specialized _track! for RungeKutta that bypasses the unpacking system.
Gets field function directly from Beamlines.field_calc and passes the full element.
"""
function _track!(
  coords::Coords,
  bunch::Bunch,
  ele::LineElement,
  tm::RungeKutta,
  ramp_without_rf;
  kwargs...
)
  # Get basic element properties
  L = float(ele.L)
  ap = deval(ele.AlignmentParams)
  bp = deval(ele.BendParams)
  dp = deval(ele.ApertureParams)
  patch = deval(ele.PatchParams)
  R_ref = bunch.R_ref

  # Setup reference state
  beta_gamma_ref = R_to_beta_gamma(bunch.species, bunch.R_ref)
  kc = KernelChain(Val{6}(), RefState(bunch.t_ref, beta_gamma_ref))

  # Evolve time through whole element
  bunch.t_ref += L/beta_gamma_to_v(beta_gamma_ref)

  # Handle reference momentum ramping
  if R_ref isa TimeDependentParam
    R_ref_initial = bunch.R_ref
    R_ref_final = R_ref(bunch.t_ref)
    if !(R_ref_initial ≈ R_ref_final)
      kc = push(kc, KernelCall(BeamTracking.update_P0!, (R_ref_initial, R_ref_final, ramp_without_rf)))
      setfield!(bunch, :R_ref, R_ref_final)
    end
  end

  if isactive(patch)
    error("RungeKutta tracking does not support patch elements")
  end

  # Entrance aperture and alignment
  if isactive(ap)
    if isactive(dp)
      if dp.aperture_shifts_with_body
        kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, true)))
        kc = push(kc, @inline(aperture(tm, bunch, dp, true)))
      else
        kc = push(kc, @inline(aperture(tm, bunch, dp, true)))
        kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, true)))
      end
    else
      kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, true)))
    end
  elseif isactive(dp)
    kc = push(kc, @inline(aperture(tm, bunch, dp, true)))
  end

    # Only track through body if element has length
    if L <= 0.0
      error("RungeKutta tracking does not support zero-length elements") 
    end

    # Setup physics parameters
    species, R_ref = bunch.species, bunch.R_ref
    tilde_m, gamsqr_0, beta_0 = BeamTracking.drift_params(species, R_ref)
    charge = chargeof(species) / BeamTracking.E_CHARGE
    p0c = BeamTracking.R_to_pc(species, R_ref)
    mc2 = massof(species)

    # Determine step size to use
    if tm.ds_step > 0
      ds_step = tm.ds_step
    elseif tm.n_steps > 0
      ds_step = L / tm.n_steps
    else
      ds_step = BeamTracking.DEFAULT_RK4_DS_STEP
    end

    s_span = (0.0, L)

    # Get curvature from BendParams if present
    g_bend = isactive(bp) ? bp.g : 0.0

    # Get field function from Beamlines and pass full element
    field_func = Beamlines.em_field_calc(ele)

    params = (beta_0, gamsqr_0, tilde_m, charge, p0c, mc2, s_span, ds_step, g_bend,
              field_func, ele)
    kc = push(kc, KernelCall(BeamTracking.RungeKuttaTracking.rk4_kernel!, params))

  # Exit aperture and alignment
  if isactive(ap)
    if isactive(dp)
      if dp.aperture_shifts_with_body
        kc = push(kc, @inline(aperture(tm, bunch, dp, false)))
        kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, false)))
      else
        kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, false)))
        kc = push(kc, @inline(aperture(tm, bunch, dp, false)))
      end
    else
      kc = push(kc, @inline(alignment(tm, bunch, ap, bp, L, false)))
    end
  elseif isactive(dp)
    kc = push(kc, @inline(aperture(tm, bunch, dp, false)))
  end

  # Launch kernels
  @noinline launch!(coords, kc; kwargs...)
  return nothing
end

# =========== ALIGNMENT AND APERTURE ============= #

@inline alignment(tm::RungeKutta, bunch, alignmentparams, bendparams, L, entering) =
  alignment(Exact(), bunch, alignmentparams, bendparams, L, entering)

@inline aperture(tm::RungeKutta, bunch, apertureparams, entering) =
  aperture(Exact(), bunch, apertureparams, entering)
