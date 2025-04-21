
@testset "LinearTracking" begin
  @testset "Utility functions" begin
    # Thick Quadrupole
    mf(K1,L) = [cos(sqrt(K1)*L)            sincu(sqrt(K1)*L)*L;  
                -sqrt(K1)*sin(sqrt(K1)*L)  cos(sqrt(K1)*L)     ]
    md(K1,L) = [cosh(sqrt(K1)*L)           sinhcu(sqrt(K1)*L)*L; 
                sqrt(K1)*sinh(sqrt(K1)*L)  cosh(sqrt(K1)*L)     ]

    L = 1.2
    K1 = 0.36
  
    # Focusing
    @test all(LinearTracking.linear_quad_matrices(K1, L) .== (mf(K1,L), md(K1,L)))

    # Defocusing
    @test all(LinearTracking.linear_quad_matrices(-K1, L) .== (md(K1,L), mf(K1,L)))


    # Thin Quadrupole
    mft(K1L) = [one(K1L) zero(K1L);
                -K1L     one(K1L)  ]
    mdt(K1L) = [one(K1L) zero(K1L);
                K1L      one(K1L)  ]
    K1L = K1*L

    @test all(LinearTracking.linear_thin_quad_matrices(K1L) .== (mft(K1L), mdt(K1L)))
    @test all(LinearTracking.linear_thin_quad_matrices(-K1L) .== (mdt(K1L), mft(K1L)))
    
  end

  @testset "Kernels" begin
    Ls = rand(Float64)
    Lt = TPS{D}(rand(Float64))

    gamma_0s = rand(Float64)
    gamma_0t = TPS{D}(rand(Float64))

    # Drift =======================================================================
    M_drift(L, gamma_0) = [1.0  L    0.0  0.0  0.0  0.0;
                           0.0  1.0  0.0  0.0  0.0  0.0;
                           0.0  0.0  1.0  L    0.0  0.0;
                           0.0  0.0  0.0  1.0  0.0  0.0;
                           0.0  0.0  0.0  0.0  1.0  L/gamma_0^2;
                           0.0  0.0  0.0  0.0  0.0  1.0] 
    # Scalar parameters
    test_matrix(LinearTracking.linear_drift!, M_drift(Ls,gamma_0s), Ls, Ls/gamma_0s^2)

    # GTPSA parameters
    test_matrix(LinearTracking.linear_drift!, M_drift(Lt,gamma_0t), Lt, Lt/gamma_0t^2)

    # Coast uncoupled  =============================================================
    function coast_uncoupled(::Type{T}) where {T}
      mx = T[1 2; 3 4]
      my = T[5 6; 7 8]
      r56 = T(9)
      d = T[10 11 12 13 14]
      t = T[15 16 17 18 19]
      return mx, my, r56, d, t
    end

    # Scalar parameters
    mxs, mys, r56s, ds, ts = coast_uncoupled(Float64)
    Ms = zeros(6,6)
    Ms[1:2, 1:2] = mxs
    Ms[3:4, 3:4] = mys
    Ms[5, 5] = 1.0
    Ms[6, 6] = 1.0
    Ms[5, 6] = r56s
    test_matrix(LinearTracking.linear_coast_uncoupled!, Ms, mxs, mys, r56s)
    Ms[5, 1:4] = ts
    Ms[6, 1:4] = ds
    test_matrix(LinearTracking.linear_coast_uncoupled!, Ms, mxs, mys, r56s, ds, ts)

    # GTPSA parameters
    mxt, myt, r56t, dt, tt = coast_uncoupled(TPS64{D})
    Mt = zeros(6,6)
    Mt[1:2, 1:2] = mxt
    Mt[3:4, 3:4] = myt
    Mt[5, 5] = 1.0
    Mt[6, 6] = 1.0
    Mt[5, 6] = r56t
    test_matrix(LinearTracking.linear_coast_uncoupled!, Mt, mxt, myt, r56t)
    Mt[5, 1:4] = tt
    Mt[6, 1:4] = dt
    test_matrix(LinearTracking.linear_coast_uncoupled!, Mt, mxt, myt, r56t, dt, tt)
  end
end
