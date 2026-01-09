@testset "RungeKuttaTracking" begin
  using BeamTracking
  using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                      RungeKuttaTracking, Bunch, STATE_ALIVE, STATE_LOST_PZ, E_CHARGE

  # Helper function to setup tracking parameters
  function setup_particle(kinetic_energy=5e3)  # 5 keV default
    species = Species("electron")
    mc2 = massof(species)  # eV
    ek = kinetic_energy
    βγ = sqrt(ek / mc2 * (ek / mc2 + 2))
    pc = mc2 * βγ
    R_ref = pc_to_R(species, pc)

    # Calculate tracking parameters
    beta_gamma_0 = R_to_beta_gamma(species, R_ref)
    tilde_m = 1 / beta_gamma_0
    gamsqr_0 = 1 + beta_gamma_0^2
    beta_0 = beta_gamma_0 / sqrt(gamsqr_0)
    charge = chargeof(species)
    p0c = R_to_pc(species, R_ref)

    return species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2
  end

  # Field functions with new signature
  function drift(x, px, y, py, z, pz, s, params)
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  end

  function uniform_efield(x, px, y, py, z, pz, s, params)
    Ex = params.Ex  # V/m
    return (Ex, 0.0, 0.0, 0.0, 0.0, 0.0)
  end

  function uniform_bfield(x, px, y, py, z, pz, s, params)
    Bz = params.Bz  # Tesla
    return (0.0, 0.0, 0.0, 0.0, 0.0, Bz)
  end

  @testset "Pure drift" begin
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    # Create bunch with small transverse momentum
    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch.coords.v[1, 1] = 0.0   # x0 = 0
    bunch.coords.v[1, 2] = 0.01  # px0 = 0.01
    bunch.coords.v[1, 3] = 0.0   # y0 = 0
    bunch.coords.v[1, 4] = 0.0   # py0 = 0
    bunch.coords.v[1, 5] = 0.0   # z0 = 0
    bunch.coords.v[1, 6] = 0.0   # pz0 = 0

    s_span = (0.0, 1.0)  # 1 meter arc length
    field_params = nothing
    ds_step = 0.01  # 1 cm step size
    g_bend = 0.0

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   drift, field_params)

    # For drift, dx/ds ≈ px (for small px and pz ≈ 0)
    # So x_final ≈ x0 + px * L
    @test isapprox(bunch.coords.v[1, 1], 0.01, rtol=1e-3)  # x ≈ 0.01 m
    @test isapprox(bunch.coords.v[1, 2], 0.01, rtol=1e-5)  # px unchanged
    @test bunch.coords.v[1, 3] ≈ 0.0  # y unchanged
    @test bunch.coords.v[1, 4] ≈ 0.0  # py unchanged
  end

  @testset "Uniform E-field - weak field" begin
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    # Start with small initial momentum in x (can't integrate arc length from rest)
    bunch.coords.v[1, 1] = 0.0
    bunch.coords.v[1, 2] = 0.001  # Small initial px
    bunch.coords.v[1, 3] = 0.0
    bunch.coords.v[1, 4] = 0.0
    bunch.coords.v[1, 5] = 0.0
    bunch.coords.v[1, 6] = 0.0

    px_initial = bunch.coords.v[1, 2]
    x_initial = bunch.coords.v[1, 1]

    s_span = (0.0, 1.0)  # 1 meter arc length
    field_params = (Ex=-1e4,)  # -10 kV/m (negative field accelerates electron in +x)
    ds_step = 0.01  # 1 cm step size
    g_bend = 0.0

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend, 
                                   uniform_efield, field_params)

    # Electron in negative E-field should accelerate in +x direction
    @test bunch.coords.v[1, 2] > px_initial  # px should increase
    @test bunch.coords.v[1, 1] > x_initial  # x should increase
    @test bunch.coords.v[1, 3] ≈ 0.0  # y unchanged
    @test bunch.coords.v[1, 4] ≈ 0.0  # py unchanged
  end

  @testset "Uniform B-field - circular motion" begin
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    # Initial velocity in x-direction
    bunch.coords.v[1, 1] = 0.0   # x0 = 0
    bunch.coords.v[1, 2] = 0.01  # px0 = 0.01
    bunch.coords.v[1, 3] = 0.0   # y0 = 0
    bunch.coords.v[1, 4] = 0.0   # py0 = 0
    bunch.coords.v[1, 5] = 0.0   # z0 = 0
    bunch.coords.v[1, 6] = 0.0   # pz0 = 0

    s_span = (0.0, 1.0)  # 1 meter
    field_params = (Bz=0.01,)  # 0.01 Tesla
    ds_step = 0.001  # 1 mm step size
    g_bend = 0.0

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   uniform_bfield, field_params)

    # In uniform B-field, particle should follow circular path
    # Total transverse momentum should be conserved
    pt2 = bunch.coords.v[1, 2]^2 + bunch.coords.v[1, 4]^2
    @test isapprox(pt2, 0.01^2, rtol=1e-4)
  end

  @testset "Particle loss detection" begin
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    # Set unphysical initial momenta (vt² > 1)
    bunch.coords.v[1, 1] = 0.0
    bunch.coords.v[1, 2] = 1.5  # px too large
    bunch.coords.v[1, 3] = 0.0
    bunch.coords.v[1, 4] = 0.0
    bunch.coords.v[1, 5] = 0.0
    bunch.coords.v[1, 6] = 0.0  # pz = 0, so rel_p = 1

    s_span = (0.0, 1.0)
    field_params = nothing
    ds_step = 0.1  # 10 cm step size
    g_bend = 0.0

    RungeKuttaTracking.rk4_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, ds_step, g_bend,
                                   drift, field_params)

    # Particle should be marked as lost
    @test bunch.coords.state[1] == STATE_LOST_PZ
  end

  @testset "Convergence test" begin
    species, R_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle()

    bunch1 = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch2 = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch1.coords.v[1, 2] = 0.01  # px0 = 0.01
    bunch2.coords.v[1, 2] = 0.01

    s_span = (0.0, 1.0)
    field_params = (Ex=1e4,)
    g_bend = 0.0

    # Track with different step sizes
    RungeKuttaTracking.rk4_kernel!(1, bunch1.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, 0.02, g_bend,
                                   uniform_efield, field_params)
    RungeKuttaTracking.rk4_kernel!(1, bunch2.coords, beta_0, gamsqr_0, tilde_m,
                                   charge, p0c, mc2, s_span, 0.005, g_bend,
                                   uniform_efield, field_params)

    # Results should be similar with finer steps being more accurate
    @test isapprox(bunch1.coords.v[1, 1], bunch2.coords.v[1, 1], rtol=1e-2)
    @test isapprox(bunch1.coords.v[1, 2], bunch2.coords.v[1, 2], rtol=1e-2)
  end

  @testset "Integration with track! and BeamlinesExt" begin
    using Beamlines

    species = Species("electron")
    mc2 = massof(species)
    ek = 5e3  # 5 keV
    βγ = sqrt(ek / mc2 * (ek / mc2 + 2))
    pc = mc2 * βγ
    R_ref = pc_to_R(species, pc)

    # Create a simple drift element
    L_drift = 1.0  # 1 meter
    drift_ele = Drift(L=L_drift)
    drift_ele.tracking_method = RungeKutta()

    # Create bunch with small transverse momentum
    bunch = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch.coords.v[1, 1] = 0.0   # x0 = 0
    bunch.coords.v[1, 2] = 0.01  # px0 = 0.01
    bunch.coords.v[1, 3] = 0.0   # y0 = 0
    bunch.coords.v[1, 4] = 0.0   # py0 = 0
    bunch.coords.v[1, 5] = 0.0   # z0 = 0
    bunch.coords.v[1, 6] = 0.0   # pz0 = 0

    # Track through drift using track!
    track!(bunch, drift_ele)

    # For drift, dx/ds ≈ px (for small px and pz ≈ 0)
    # So x_final ≈ x0 + px * L
    @test isapprox(bunch.coords.v[1, 1], 0.01, rtol=1e-3)  # x ≈ 0.01 m
    @test isapprox(bunch.coords.v[1, 2], 0.01, rtol=1e-5)  # px unchanged
    @test bunch.coords.v[1, 3] ≈ 0.0  # y unchanged
    @test bunch.coords.v[1, 4] ≈ 0.0  # py unchanged
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "RungeKutta with different step configurations" begin
    using Beamlines

    species = Species("electron")
    mc2 = massof(species)
    ek = 5e3
    βγ = sqrt(ek / mc2 * (ek / mc2 + 2))
    pc = mc2 * βγ
    R_ref = pc_to_R(species, pc)

    L_drift = 1.0

    # Test with ds_step
    drift_ds = Drift(L=L_drift)
    drift_ds.tracking_method = RungeKutta(ds_step=0.1)
    bunch_ds = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch_ds.coords.v[1, 2] = 0.01
    track!(bunch_ds, drift_ds)

    # Test with n_steps
    drift_ns = Drift(L=L_drift)
    drift_ns.tracking_method = RungeKutta(ds_step=-1.0, n_steps=50)
    bunch_ns = Bunch(zeros(1, 6), R_ref=R_ref, species=species)
    bunch_ns.coords.v[1, 2] = 0.01
    track!(bunch_ns, drift_ns)

    # Both should give similar results
    @test isapprox(bunch_ds.coords.v[1, 1], bunch_ns.coords.v[1, 1], rtol=1e-2)
    @test isapprox(bunch_ds.coords.v[1, 2], bunch_ns.coords.v[1, 2], rtol=1e-4)
  end

end
