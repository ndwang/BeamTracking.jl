@testset "Beamlines" begin
  include("lattices/esr.jl")


  @testset "Linear" begin
    d = Descriptor(6, 1)
    b0 = Bunch(collect(transpose(@vars(d))), Brho_ref=ring.Brho_ref)
    foreach(t->t.tracking_method=Linear(), ring.line)
    track!(b0, ring)
    M_ESR = [  0.8763088913632391E+00  0.2842332738570903E+00 -0.9233408598814828E-06 -0.1104742931103878E-06  0.0000000000000000E+00 -0.8843595261589663E-07
              -0.8165279324836310E+00  0.8763069736287854E+00 -0.1898521129265218E-05 -0.1113630193834745E-06  0.0000000000000000E+00 -0.1417461685411299E-06
               0.1352460571822227E-07 -0.2969583000777938E-07  0.6374265424413608E+00  0.4391919124687460E-01  0.0000000000000000E+00  0.5592919170820314E-14
              -0.4079601917518316E-05  0.1052556482124766E-05 -0.1351773220872402E+02  0.6374258159142916E+00  0.0000000000000000E+00 -0.4735132312703365E-13
               0.1964238541540386E-06 -0.3720806396318763E-07 -0.8403098936929155E-14 -0.1661103086541925E-15  0.1000000000000000E+01 -0.2358669003370585E+01
               0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00  0.1000000000000000E+01  ]

    @test GTPSA.jacobian(b0.v) ≈ M_ESR

    bblring = BitsBeamline(ring, store_normalized=true)
    b0 = Bunch(collect(transpose(@vars(d))), Brho_ref=ring.Brho_ref)
    track!(b0, bblring)
    @test GTPSA.jacobian(b0.v) ≈ M_ESR
  end

  @testset "Exact" begin
    # Define elements
    ele_patch = LineElement(dt=1e-9, dx=2.0, dy=3.0, dz=4.0, dx_rot=-5.0, dy_rot=6.0, dz_rot=7.0, L=-1.9458360380198412, tracking_method=Exact())
    ele_drift = LineElement(L=1.0, tracking_method=Exact())
    ele_sol = LineElement(L=1.0, Ks=2.0, tracking_method=Exact())
    ele_kick = LineElement(L=1.0, K0L=1.0, tracking_method=Exact())
    ele_bend = LineElement(L=1.0, g=0.01, tracking_method=Exact())
    ele_patch_bend = LineElement(L=1.0, g=0.01, dy=3.0, dz_rot=0.3, tracking_method=Exact())
    ele_patch_sol = LineElement(L=1.0, Ks=1.0, dt=1.0, tracking_method=Exact())
    ele_bend_quad = LineElement(L=1.0, g=0.01, K1=1.0, tracking_method=Exact())

    v = zeros(1,6)
    bunch = Bunch(v, species=ELECTRON, Brho_ref=1.0)

    @testset "patch only" begin
        test_map(track!, "bmad_maps/patch.jl", kernel_test=false, ele=ele_patch, p0c=10e6, species=ELECTRON, tol=5e-10)
    end
    @testset "drift" begin
      test_map(track!, "bmad_maps/drift.jl", kernel_test=false, ele=ele_drift, p0c=10e6, species=ELECTRON, tol=5e-10)
    end
    @testset "thick solenoid" begin
      test_map(track!, "bmad_maps/solenoid.jl", kernel_test=false, ele=ele_sol, p0c=10e6, species=ELECTRON, tol=5e-10)
    end
    @testset "kick" begin
      @test_throws ErrorException track!(bunch, ele_kick)
    end
    @testset "bend only" begin
        @test_throws ErrorException track!(bunch, ele_bend)
    end
    @testset "patch with bend" begin
        @test_throws ErrorException track!(bunch, ele_patch_bend)
    end
    @testset "patch with multipole" begin
        @test_throws ErrorException track!(bunch, ele_patch_sol)
    end
    @testset "combined function bend" begin
        @test_throws ErrorException track!(bunch, ele_bend_quad)
    end
  end
end