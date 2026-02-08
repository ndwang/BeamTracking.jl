@testset "FieldMap" begin
  using BeamTracking
  using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                      RungeKuttaTracking, Bunch, STATE_ALIVE, STATE_LOST_PZ, E_CHARGE, C_LIGHT,
                      RectGrid3D, CylGrid2D, FieldMap,
                      interp_component, interp_field, fieldmap_em_field,
                      grid_min
  using StaticArrays

  # Helper to set up tracking parameters
  function setup_particle(pc=1e9)
    species = Species("electron")
    mc2 = massof(species)
    p_over_q_ref = pc_to_R(species, pc)
    beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
    tilde_m = 1 / beta_gamma_0
    gamsqr_0 = 1 + beta_gamma_0^2
    beta_0 = beta_gamma_0 / sqrt(gamsqr_0)
    charge = chargeof(species)
    p0c = R_to_pc(species, p_over_q_ref)
    return species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2
  end

  # ================================================================
  # Phase 1: Grid Types, FieldMap Struct, Interpolation
  # ================================================================

  @testset "Grid types" begin
    g3 = RectGrid3D((0.0, 0.0, 0.0), (0.1, 0.2, 0.3), (0, 0, 0), (5, 6, 7))
    @test grid_min(g3, 1) == 0.0
    @test grid_min(g3, 2) == 0.0
    @test grid_min(g3, 3) == 0.0

    # With offset and lower bound
    g3b = RectGrid3D((-1.0, -2.0, 0.5), (0.1, 0.2, 0.3), (1, 2, 0), (5, 6, 7))
    @test grid_min(g3b, 1) == -1.0 + 1 * 0.1
    @test grid_min(g3b, 2) == -2.0 + 2 * 0.2
    @test grid_min(g3b, 3) == 0.5

    g2 = CylGrid2D((0.0, 0.0), (0.01, 0.1), (0, 0), (10, 20))
    @test grid_min(g2, 1) == 0.0
    @test grid_min(g2, 2) == 0.0
  end

  @testset "FieldMap constructor - RectGrid3D" begin
    grid = RectGrid3D((0.0, 0.0, 0.0), (0.1, 0.1, 0.1), (0, 0, 0), (3, 3, 3))

    # B-only map
    B = ones(3, 3, 3, 3)
    fm = FieldMap(grid; B=B)
    @test fm.B_re === B
    @test fm.E_re === nothing
    @test fm.E_im === nothing
    @test fm.B_im === nothing
    @test fm.frequency == 0.0
    @test fm.harmonic == 0
    @test fm.eleAnchorPt == :beginning

    # E-only map
    E = ones(3, 3, 3, 3)
    fm_e = FieldMap(grid; E=E)
    @test fm_e.E_re === E
    @test fm_e.B_re === nothing

    # Both E and B
    fm_both = FieldMap(grid; E=E, B=B)
    @test fm_both.E_re === E
    @test fm_both.B_re === B

    # Wrong shape should error
    B_bad = ones(3, 4, 3, 3)
    @test_throws AssertionError FieldMap(grid; B=B_bad)

    # Neither E nor B should error
    @test_throws ErrorException FieldMap(grid)
  end

  @testset "FieldMap constructor - CylGrid2D" begin
    grid = CylGrid2D((0.0, 0.0), (0.01, 0.1), (0, 0), (5, 10))

    B = ones(3, 5, 10)
    fm = FieldMap(grid; B=B)
    @test fm.B_re === B
    @test fm.E_re === nothing

    # Wrong shape
    B_bad = ones(3, 4, 10)
    @test_throws AssertionError FieldMap(grid; B=B_bad)

    # Neither
    @test_throws ErrorException FieldMap(grid)
  end

  @testset "Trilinear interpolation - accuracy" begin
    # Create a grid with known analytical function: f(x,y,z) = 2x + 3y + 5z
    # Linear function should be interpolated exactly by trilinear
    nx, ny, nz = 5, 6, 7
    dx, dy, dz = 0.1, 0.2, 0.3
    grid = RectGrid3D((0.0, 0.0, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    arr = zeros(3, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
      x = (ix - 1) * dx
      y = (iy - 1) * dy
      z = (iz - 1) * dz
      arr[1, ix, iy, iz] = 2x + 3y + 5z        # Component 1
      arr[2, ix, iy, iz] = x^2                   # Component 2 (nonlinear for error test)
      arr[3, ix, iy, iz] = 7.0                   # Component 3 (constant)
    end

    fm = FieldMap(grid; B=arr)

    # Test at off-grid point (linear function should be exact)
    x_test, y_test, z_test = 0.15, 0.35, 0.45
    val1 = interp_component(1, arr, grid, x_test, y_test, z_test)
    @test val1 ≈ 2*x_test + 3*y_test + 5*z_test atol=1e-14

    # Constant function should be exact
    val3 = interp_component(3, arr, grid, x_test, y_test, z_test)
    @test val3 ≈ 7.0 atol=1e-14

    # Vector interpolation
    f1, f2, f3 = interp_field(arr, grid, x_test, y_test, z_test)
    @test f1 ≈ 2*x_test + 3*y_test + 5*z_test atol=1e-14
    @test f3 ≈ 7.0 atol=1e-14
  end

  @testset "Trilinear interpolation - on-grid-point" begin
    nx, ny, nz = 4, 4, 4
    grid = RectGrid3D((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0, 0, 0), (nx, ny, nz))

    arr = randn(3, nx, ny, nz)

    # Test at every grid point
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
      x = Float64(ix - 1)
      y = Float64(iy - 1)
      z = Float64(iz - 1)
      for c in 1:3
        @test interp_component(c, arr, grid, x, y, z) ≈ arr[c, ix, iy, iz] atol=1e-14
      end
    end
  end

  @testset "Trilinear interpolation - out-of-bounds" begin
    nx, ny, nz = 3, 3, 3
    grid = RectGrid3D((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0, 0, 0), (nx, ny, nz))
    arr = ones(3, nx, ny, nz)

    # Test all 6 faces: outside should be zero
    @test interp_component(1, arr, grid, -0.1, 0.5, 0.5) == 0.0
    @test interp_component(1, arr, grid,  2.1, 0.5, 0.5) == 0.0
    @test interp_component(1, arr, grid,  0.5, -0.1, 0.5) == 0.0
    @test interp_component(1, arr, grid,  0.5,  2.1, 0.5) == 0.0
    @test interp_component(1, arr, grid,  0.5, 0.5, -0.1) == 0.0
    @test interp_component(1, arr, grid,  0.5, 0.5,  2.1) == 0.0

    # Inside should not be zero
    @test interp_component(1, arr, grid, 0.5, 0.5, 0.5) ≈ 1.0
  end

  @testset "Trilinear interpolation - single-cell grid" begin
    grid = RectGrid3D((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0, 0, 0), (2, 2, 2))
    arr = zeros(3, 2, 2, 2)
    arr[1, 1, 1, 1] = 1.0
    arr[1, 2, 2, 2] = 1.0

    @test interp_component(1, arr, grid, 0.0, 0.0, 0.0) ≈ 1.0
    @test interp_component(1, arr, grid, 1.0, 1.0, 1.0) ≈ 1.0
    # Center: only corner [1,1,1] and [2,2,2] contribute with weight 1/8 each
    @test interp_component(1, arr, grid, 0.5, 0.5, 0.5) ≈ 0.25
  end

  @testset "Bilinear interpolation - accuracy" begin
    # Linear function f(r,z) = 2r + 5z — should be exact
    nr, nz = 5, 7
    dr, dz = 0.01, 0.1
    grid = CylGrid2D((0.0, 0.0), (dr, dz), (0, 0), (nr, nz))

    arr = zeros(3, nr, nz)
    for ir in 1:nr, iz in 1:nz
      r = (ir - 1) * dr
      z = (iz - 1) * dz
      arr[1, ir, iz] = 2r + 5z
      arr[2, ir, iz] = r^2
      arr[3, ir, iz] = 3.0
    end

    r_test, z_test = 0.015, 0.35
    val1 = interp_component(1, arr, grid, r_test, z_test)
    @test val1 ≈ 2*r_test + 5*z_test atol=1e-14

    val3 = interp_component(3, arr, grid, r_test, z_test)
    @test val3 ≈ 3.0 atol=1e-14

    f1, f2, f3 = interp_field(arr, grid, r_test, z_test)
    @test f1 ≈ 2*r_test + 5*z_test atol=1e-14
    @test f3 ≈ 3.0 atol=1e-14
  end

  @testset "Bilinear interpolation - on-grid-point" begin
    nr, nz = 4, 5
    grid = CylGrid2D((0.0, 0.0), (1.0, 1.0), (0, 0), (nr, nz))
    arr = randn(3, nr, nz)

    for ir in 1:nr, iz in 1:nz
      r = Float64(ir - 1)
      z = Float64(iz - 1)
      for c in 1:3
        @test interp_component(c, arr, grid, r, z) ≈ arr[c, ir, iz] atol=1e-14
      end
    end
  end

  @testset "Bilinear interpolation - out-of-bounds" begin
    nr, nz = 3, 3
    grid = CylGrid2D((0.0, 0.0), (1.0, 1.0), (0, 0), (nr, nz))
    arr = ones(3, nr, nz)

    @test interp_component(1, arr, grid, -0.1, 0.5) == 0.0
    @test interp_component(1, arr, grid,  2.1, 0.5) == 0.0
    @test interp_component(1, arr, grid,  0.5, -0.1) == 0.0
    @test interp_component(1, arr, grid,  0.5,  2.1) == 0.0

    @test interp_component(1, arr, grid, 0.5, 0.5) ≈ 1.0
  end

  @testset "Type stability - interpolation" begin
    grid3 = RectGrid3D((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0, 0, 0), (3, 3, 3))
    arr3 = ones(3, 3, 3, 3)
    @inferred interp_component(1, arr3, grid3, 0.5, 0.5, 0.5)
    @inferred interp_field(arr3, grid3, 0.5, 0.5, 0.5)

    grid2 = CylGrid2D((0.0, 0.0), (1.0, 1.0), (0, 0), (3, 3))
    arr2 = ones(3, 3, 3)
    @inferred interp_component(1, arr2, grid2, 0.5, 0.5)
    @inferred interp_field(arr2, grid2, 0.5, 0.5)
  end

  # ================================================================
  # Phase 2: Field Evaluation
  # ================================================================

  @testset "Rectangular field evaluation - B only" begin
    nx, ny, nz = 5, 5, 10
    dx, dy, dz = 0.01, 0.01, 0.1
    grid = RectGrid3D((0.0, 0.0, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    # Uniform Bz = 1.0 T
    B = zeros(3, nx, ny, nz)
    B[3, :, :, :] .= 1.0
    fm = FieldMap(grid; B=B)

    # E should be zero, Bz should be 1.0
    Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(0.02, 0.02, 0.0, 0.0, 0.5, fm, 0.0)
    @test Ex == 0.0
    @test Ey == 0.0
    @test Ez == 0.0
    @test Bx ≈ 0.0 atol=1e-15
    @test By ≈ 0.0 atol=1e-15
    @test Bz ≈ 1.0 atol=1e-14
  end

  @testset "Rectangular field evaluation - E only" begin
    nx, ny, nz = 5, 5, 10
    dx, dy, dz = 0.01, 0.01, 0.1
    grid = RectGrid3D((0.0, 0.0, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    # Uniform Ez = 1e6 V/m
    E = zeros(3, nx, ny, nz)
    E[3, :, :, :] .= 1e6
    fm = FieldMap(grid; E=E)

    Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(0.02, 0.02, 0.0, 0.0, 0.5, fm, 0.0)
    @test Bx == 0.0
    @test By == 0.0
    @test Bz == 0.0
    @test Ex ≈ 0.0 atol=1e-10
    @test Ey ≈ 0.0 atol=1e-10
    @test Ez ≈ 1e6 atol=1e-8
  end

  @testset "Rectangular field evaluation - z_offset" begin
    nx, ny, nz = 3, 3, 5
    grid = RectGrid3D((0.0, 0.0, 0.0), (0.01, 0.01, 0.25), (0, 0, 0), (nx, ny, nz))

    # B varies linearly with z: Bz(z) = z
    B = zeros(3, nx, ny, nz)
    for iz in 1:nz
      B[3, :, :, iz] .= (iz - 1) * 0.25
    end
    fm = FieldMap(grid; B=B)

    # At s=0.5 with z_offset=0 -> z_map=0.5, expect Bz=0.5
    _, _, _, _, _, Bz = fieldmap_em_field(0.0, 0.0, 0.0, 0.0, 0.5, fm, 0.0)
    @test Bz ≈ 0.5 atol=1e-14

    # At s=0.5 with z_offset=-0.25 -> z_map=0.25
    _, _, _, _, _, Bz2 = fieldmap_em_field(0.0, 0.0, 0.0, 0.0, 0.5, fm, -0.25)
    @test Bz2 ≈ 0.25 atol=1e-14
  end

  @testset "Cylindrical field evaluation - solenoid" begin
    # Uniform Bz solenoid: Br=0, Btheta=0, Bz=const
    nr, nz_grid = 5, 10
    grid = CylGrid2D((0.0, 0.0), (0.01, 0.1), (0, 0), (nr, nz_grid))

    B = zeros(3, nr, nz_grid)
    B[3, :, :] .= 1.0  # Bz = 1.0 T
    fm = FieldMap(grid; B=B)

    # At (x=0.02, y=0.01, z=0.0)
    Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(0.02, 0.01, 0.0, 0.0, 0.5, fm, 0.0)
    @test Ex == 0.0
    @test Ey == 0.0
    @test Ez == 0.0
    @test Bx ≈ 0.0 atol=1e-15
    @test By ≈ 0.0 atol=1e-15
    @test Bz ≈ 1.0 atol=1e-14
  end

  @testset "Cylindrical field evaluation - r=0 edge case" begin
    nr, nz_grid = 5, 5
    grid = CylGrid2D((0.0, 0.0), (0.01, 0.1), (0, 0), (nr, nz_grid))

    B = zeros(3, nr, nz_grid)
    B[3, :, :] .= 1.0
    fm = FieldMap(grid; B=B)

    # At r=0 exactly: should not produce NaN/Inf
    Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(0.0, 0.0, 0.0, 0.0, 0.2, fm, 0.0)
    @test isfinite(Ex)
    @test isfinite(Ey)
    @test isfinite(Ez)
    @test isfinite(Bx)
    @test isfinite(By)
    @test isfinite(Bz)
    @test Bz ≈ 1.0 atol=1e-14
  end

  @testset "Cylindrical field evaluation - radial field" begin
    # Pure radial field: Fr = r, Ftheta = 0, Fz = 0
    # In Cartesian: Fx = Fr*cos(theta) = x, Fy = Fr*sin(theta) = y
    nr, nz_grid = 10, 5
    dr = 0.01
    grid = CylGrid2D((0.0, 0.0), (dr, 0.1), (0, 0), (nr, nz_grid))

    B = zeros(3, nr, nz_grid)
    for ir in 1:nr
      r = (ir - 1) * dr
      B[1, ir, :] .= r  # Br = r
    end
    fm = FieldMap(grid; B=B)

    x_test, y_test = 0.03, 0.04  # r = 0.05
    Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(x_test, y_test, 0.0, 0.0, 0.2, fm, 0.0)
    # Bx = Br*cos(theta) = r * x/r = x
    @test Bx ≈ x_test atol=1e-10
    # By = Br*sin(theta) = r * y/r = y
    @test By ≈ y_test atol=1e-10
    @test Bz ≈ 0.0 atol=1e-15
  end

  @testset "Type stability - field evaluation" begin
    # Rectangular B-only
    grid3 = RectGrid3D((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0, 0, 0), (3, 3, 3))
    B3 = ones(3, 3, 3, 3)
    fm3_b = FieldMap(grid3; B=B3)
    @inferred fieldmap_em_field(0.5, 0.5, 0.0, 0.0, 0.5, fm3_b, 0.0)

    # Rectangular E-only
    E3 = ones(3, 3, 3, 3)
    fm3_e = FieldMap(grid3; E=E3)
    @inferred fieldmap_em_field(0.5, 0.5, 0.0, 0.0, 0.5, fm3_e, 0.0)

    # Rectangular E+B
    fm3_eb = FieldMap(grid3; E=E3, B=B3)
    @inferred fieldmap_em_field(0.5, 0.5, 0.0, 0.0, 0.5, fm3_eb, 0.0)

    # Cylindrical B-only
    grid2 = CylGrid2D((0.0, 0.0), (1.0, 1.0), (0, 0), (3, 3))
    B2 = ones(3, 3, 3)
    fm2_b = FieldMap(grid2; B=B2)
    @inferred fieldmap_em_field(0.5, 0.5, 0.0, 0.0, 0.5, fm2_b, 0.0)
  end

  # ================================================================
  # Phase 3: RK4 Kernel Integration -- Physics Validation
  # ================================================================

  @testset "Uniform Bz field - circular orbit" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    # Create a field map with uniform Bz
    Bz_phys = 0.01  # Tesla
    nx, ny, nz = 5, 5, 20
    x_range = 0.02
    y_range = 0.02
    z_range = 1.0
    dx = 2*x_range / (nx - 1)
    dy = 2*y_range / (ny - 1)
    dz = z_range / (nz - 1)
    grid = RectGrid3D((-x_range, -y_range, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    B = zeros(3, nx, ny, nz)
    B[3, :, :, :] .= Bz_phys
    fm = FieldMap(grid; B=B)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    s_span = (0.0, z_range)
    ds_step = 0.01
    z_offset = 0.0

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, s_span, ds_step, fm, z_offset)

    # Total transverse momentum should be conserved (uniform B)
    pt2 = bunch.coords.v[1, 2]^2 + bunch.coords.v[1, 4]^2
    @test isapprox(pt2, 0.01^2, rtol=1e-4)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Uniform Bz field map vs multipole solenoid" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    Bz_phys = 0.01  # Tesla
    L = 1.0

    # --- Field map version ---
    nx, ny, nz = 5, 5, 20
    x_range = 0.02
    y_range = 0.02
    dx = 2*x_range / (nx - 1)
    dy = 2*y_range / (ny - 1)
    dz = L / (nz - 1)
    grid = RectGrid3D((-x_range, -y_range, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    B = zeros(3, nx, ny, nz)
    B[3, :, :, :] .= Bz_phys
    fm = FieldMap(grid; B=B)

    bunch_fm = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_fm.coords.v[1, BeamTracking.PXI] = 0.01

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch_fm.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.01, fm, 0.0)

    # --- Multipole version (now also uses physical fields) ---
    Bz_normalized = Bz_phys / p_over_q_ref
    mm = SVector(0)
    kn = SVector(Bz_normalized)
    ks = SVector(0.0)

    bunch_mp = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_mp.coords.v[1, BeamTracking.PXI] = 0.01

    RungeKuttaTracking.rk4_kernel!(1, bunch_mp.coords, beta_0, gamsqr_0, tilde_m,
                                    charge, p0c, mc2, (0.0, L), 0.01, 0.0, mm, kn, ks, p_over_q_ref)

    # Field map and multipole should agree closely (within interpolation error)
    @test isapprox(bunch_fm.coords.v, bunch_mp.coords.v, rtol=1e-6)
  end

  @testset "Uniform Ez field - energy gain" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    # Uniform Ez = 1e6 V/m over 1 m
    Ez_phys = 1e6  # V/m
    L = 1.0

    nx, ny, nz = 3, 3, 20
    x_range = 0.02
    y_range = 0.02
    dx = 2*x_range / (nx - 1)
    dy = 2*y_range / (ny - 1)
    dz = L / (nz - 1)
    grid = RectGrid3D((-x_range, -y_range, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    E_arr = zeros(3, nx, ny, nz)
    E_arr[3, :, :, :] .= Ez_phys
    fm = FieldMap(grid; E=E_arr)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.01, fm, 0.0)

    # Energy change for ultra-relativistic on-axis particle: dpz ~ charge * Ez * L / p0c
    dpz = bunch.coords.v[1, BeamTracking.PZI]
    expected_dpz = charge * Ez_phys * L / p0c
    @test isapprox(dpz, expected_dpz, rtol=1e-2)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Quadrupole field map vs multipole" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    # Quadrupole: Bx = G*y, By = G*x (normal quad)
    G = 10.0  # T/m gradient
    L = 0.5

    # Field map version
    nx, ny, nz = 11, 11, 20
    x_range = 0.05
    y_range = 0.05
    dx = 2*x_range / (nx - 1)
    dy = 2*y_range / (ny - 1)
    dz = L / (nz - 1)
    grid = RectGrid3D((-x_range, -y_range, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    B = zeros(3, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny
      x = -x_range + (ix - 1) * dx
      y = -y_range + (iy - 1) * dy
      B[1, ix, iy, :] .= G * y   # Bx = G*y
      B[2, ix, iy, :] .= G * x   # By = G*x
    end
    fm = FieldMap(grid; B=B)

    bunch_fm = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_fm.coords.v[1, BeamTracking.XI] = 0.005   # 5 mm offset
    bunch_fm.coords.v[1, BeamTracking.PXI] = 0.001

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch_fm.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.005, fm, 0.0)

    # Multipole version
    G_normalized = G / p_over_q_ref
    mm = SVector(2)
    kn = SVector(G_normalized)
    ks = SVector(0.0)

    bunch_mp = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch_mp.coords.v[1, BeamTracking.XI] = 0.005
    bunch_mp.coords.v[1, BeamTracking.PXI] = 0.001

    RungeKuttaTracking.rk4_kernel!(1, bunch_mp.coords, beta_0, gamsqr_0, tilde_m,
                                    charge, p0c, mc2, (0.0, L), 0.005, 0.0, mm, kn, ks, p_over_q_ref)

    # Should agree within interpolation error
    @test isapprox(bunch_fm.coords.v, bunch_mp.coords.v, rtol=1e-3)
  end

  @testset "Dead particles not updated" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    nx, ny, nz = 3, 3, 10
    grid = RectGrid3D((-0.01, -0.01, 0.0), (0.01, 0.01, 0.1), (0, 0, 0), (nx, ny, nz))
    B = zeros(3, nx, ny, nz)
    B[3, :, :, :] .= 1.0
    fm = FieldMap(grid; B=B)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 1.5  # Unphysical
    v_before = copy(bunch.coords.v)

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, 0.9), 0.1, fm, 0.0)

    @test bunch.coords.state[1] == STATE_LOST_PZ
    @test bunch.coords.v ≈ v_before
  end

  @testset "eleAnchorPt z_offset computation" begin
    L = 2.0

    grid = RectGrid3D((0.0, 0.0, 0.0), (0.01, 0.01, 0.1), (0, 0, 0), (3, 3, 25))
    B = ones(3, 3, 3, 25)

    fm_begin = FieldMap(grid; B=B, eleAnchorPt=:beginning)
    fm_center = FieldMap(grid; B=B, eleAnchorPt=:center)
    fm_end = FieldMap(grid; B=B, eleAnchorPt=:end)

    @test fm_begin.eleAnchorPt == :beginning
    @test fm_center.eleAnchorPt == :center
    @test fm_end.eleAnchorPt == :end

    # z_offset computation (done at tracking time, not in FieldMap)
    z_off_begin = 0.0
    z_off_center = -L/2
    z_off_end = -L
    @test z_off_begin == 0.0
    @test z_off_center == -1.0
    @test z_off_end == -2.0
  end

  @testset "Cylindrical solenoid - Larmor rotation via fieldmap kernel" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    Bz_phys = 0.01
    L = 1.0

    nr, nz_grid = 10, 20
    dr = 0.005
    dz = L / (nz_grid - 1)
    grid = CylGrid2D((0.0, 0.0), (dr, dz), (0, 0), (nr, nz_grid))

    B = zeros(3, nr, nz_grid)
    B[3, :, :] .= Bz_phys
    fm = FieldMap(grid; B=B)

    bunch = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch.coords.v[1, BeamTracking.PXI] = 0.01

    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.01, fm, 0.0)

    # Transverse momentum should be conserved
    pt2 = bunch.coords.v[1, 2]^2 + bunch.coords.v[1, 4]^2
    @test isapprox(pt2, 0.01^2, rtol=1e-4)
    @test bunch.coords.state[1] == STATE_ALIVE
  end

  @testset "Convergence test - fieldmap" begin
    species, p_over_q_ref, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2 = setup_particle(1e9)

    nx, ny, nz = 5, 5, 20
    x_range = 0.02
    y_range = 0.02
    L = 1.0
    dx = 2*x_range / (nx - 1)
    dy = 2*y_range / (ny - 1)
    dz = L / (nz - 1)
    grid = RectGrid3D((-x_range, -y_range, 0.0), (dx, dy, dz), (0, 0, 0), (nx, ny, nz))

    B = zeros(3, nx, ny, nz)
    B[3, :, :, :] .= 0.01
    fm = FieldMap(grid; B=B)

    # Coarse step
    bunch1 = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch1.coords.v[1, BeamTracking.PXI] = 0.01
    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch1.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.1, fm, 0.0)

    # Fine step
    bunch2 = Bunch(zeros(1, 6), p_over_q_ref=p_over_q_ref, species=species)
    bunch2.coords.v[1, BeamTracking.PXI] = 0.01
    RungeKuttaTracking.rk4_fieldmap_kernel!(1, bunch2.coords, beta_0, gamsqr_0, tilde_m,
                                             charge, p0c, mc2, (0.0, L), 0.05, fm, 0.0)

    @test isapprox(bunch1.coords.v, bunch2.coords.v, rtol=1e-2)
  end

end
