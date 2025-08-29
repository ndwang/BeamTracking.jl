using ..SpaceChargeIntegrationTracking: @sc_wrap

@inline function thick_bend_no_field(tm::SpaceChargeIntegration, bunch, bendparams, scp, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = ExactTracking.w_quaternion(0,0,-tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,-tilt)
  theta = g * L
  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, beta_0, e1, e2, theta, g, w, w_inv)
  return integration_launcher!(@sc_wrap(ExactTracking.exact_curved_drift!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function drift(tm::SpaceChargeIntegration, bunch, scp, L)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, gamsqr_0, beta_0)
  return integration_launcher!(@sc_wrap(ExactTracking.exact_drift!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_pure_bsolenoid(tm::SpaceChargeIntegration, bunch, bm0, scp, L)
  Ksol, _ = get_strengths(bm0, L, bunch.R_ref)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, gamsqr_0, beta_0, Ksol)
  return integration_launcher!(@sc_wrap(ExactTracking.exact_solenoid!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_bend_pure_bdipole(tm::SpaceChargeIntegration, bunch, bendparams, bm1, scp, L)
  g = bendparams.g_ref
  tilt = bendparams.tilt_ref
  e1 = bendparams.e1
  e2 = bendparams.e2
  w = ExactTracking.w_quaternion(0,0,-tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,-tilt)
  theta = g * L
  Kn0, Ks0 = get_strengths(bm1, L, bunch.R_ref)
  Ks0 ≈ 0 || error("A skew dipole field cannot be used in an exact bend")
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, beta_0, e1, e2, theta, g, Kn0, w, w_inv)
  return integration_launcher!(@sc_wrap(ExactTracking.exact_bend!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_bmultipole(tm::SpaceChargeIntegration, bunch, bm, scp, L)
  mm = bm.order
  kn, ks = get_strengths(bm, L, bunch.R_ref)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, gamsqr_0, beta_0, BeamTracking.anom(bunch.species), mm, kn, ks)
  return integration_launcher!(@sc_wrap(IntegrationTracking.dkd_multipole!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_bsolenoid(tm::SpaceChargeIntegration, bunch, bm, scp, L) 
  mm = bm.order
  kn, ks = get_strengths(bm, L, bunch.R_ref)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, gamsqr_0, beta_0, BeamTracking.anom(bunch.species), kn[1], mm, kn, ks)
  return integration_launcher!(@sc_wrap(IntegrationTracking.sks_multipole!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_bend_bdipole(tm::SpaceChargeIntegration, bunch, bm, scp, L)
  mm = bm.order
  kn, ks = get_strengths(bm, L, bunch.R_ref)
  k0 = sqrt(kn[1]^2 + ks[1]^2)
  tilt = atan(ks[1], kn[1])
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, beta_0, BeamTracking.anom(bunch.species), 0, 0, 0, w, w_inv, k0, mm, kn, ks)
  return integration_launcher!(@sc_wrap(IntegrationTracking.bkb_multipole!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end

@inline function thick_bquadrupole(tm::SpaceChargeIntegration, bunch, bm, scp, L)
  mm = bm.order
  kn, ks = get_strengths(bm, L, bunch.R_ref)
  k1 = sqrt(kn[1]^2 + ks[1]^2)
  tilt = atan(ks[1], kn[1]) / 2
  w = ExactTracking.w_quaternion(0,0,tilt)
  w_inv = ExactTracking.w_inv_quaternion(0,0,tilt)
  tilde_m, gamsqr_0, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
  params = (tilde_m, gamsqr_0, beta_0, BeamTracking.anom(bunch.species), w, w_inv, k1, mm, kn, ks)
  return integration_launcher!(@sc_wrap(IntegrationTracking.mkm_quadrupole!, params, beta_0, bunch.R_ref, scp, L), (), tm, L)
end