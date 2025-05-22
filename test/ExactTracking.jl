@testset "ExactTracking" begin
    @testset "Utility functions" begin
        dx_rot = 0.1
        dy_rot = -0.1
        dz_rot = 0.2
  
        W = [cos(dy_rot) 0 sin(dy_rot); 0 1 0; -sin(dy_rot) 0 cos(dy_rot)] *
          [1 0 0; 0 cos(dx_rot) sin(dx_rot); 0 -sin(dx_rot) cos(dx_rot)] *
          [cos(dz_rot) -sin(dz_rot) 0; sin(dz_rot) cos(dz_rot) 0; 0 0 1]
  
        # Test w_matrix function
        @test all(ExactTracking.w_matrix(dx_rot, dy_rot, dz_rot) .== W)
  
        Winv = [cos(dz_rot) sin(dz_rot) 0; -sin(dz_rot) cos(dz_rot) 0; 0 0 1] *
             [1 0 0; 0 cos(dx_rot) -sin(dx_rot); 0 sin(dx_rot) cos(dx_rot)] *
             [cos(dy_rot) 0 -sin(dy_rot); 0 1 0; sin(dy_rot) 0 cos(dy_rot)]
  
        # Test w_inv_matrix function
        @test all(ExactTracking.w_inv_matrix(dx_rot, dy_rot, dz_rot) .== Winv)
    end
  
    @testset "Kernels" begin
        function coeffs_approx_equal(v1, v2, ϵ)
            n = GTPSA.numcoefs(v1[1])
            for i in 1:6
                for j in 0:n-1
                    c1, c2 = v1[i][j], v2[i][j]
                    @test abs(c1 - c2) <= max(1e-12, ϵ * (abs(c1) + abs(c2)))
                end
            end
            return true
        end

        # Bmad matrices (patch, drift)
        include("bmad_maps/matrices.jl")
    
        function patch(::Type{T}) where {T}
            dx = T(1)
            dy = T(2)
            dz = T(3)
            winv = ExactTracking.w_inv_matrix(T(4),T(5),T(6))
            return dx, dy, dz, winv
        end

        # Scalar patch jacobian
        dx, dy, dz, winv = patch(Float64)
        test_matrix(ExactTracking.patch!, M_patch_bmad, 10e6, 510998.95, dx, dy, dz, winv; tol=1e-7)

        # GTPSA patch jacobian
        dx, dy, dz, winv = patch(TPS64{D})
        test_matrix(ExactTracking.patch!, M_patch_bmad, 10e6, 510998.95, dx, dy, dz, winv; tol=1e-7)

        # GTPSA patch map comparison
        include("bmad_maps/patch.jl")
        dx, dy, dz, winv = patch(TPS64{d_z})
        v = transpose(@vars(d_z))
        v = ExactTracking.patch!(1, v, zeros(eltype(v),1,9), 10e6, 510998.95,dx,dy,dz,winv)

        coeffs_approx_equal(v_z, v, 5e-8)

        # Scalar drift jacobian
        L = Float64(1)
        test_matrix(ExactTracking.exact_drift!, M_drift_bmad, L, 10e6, 510998.95; tol=1e-7)

        # GTPSA drift jacobian
        L = TPS64{D}(1)
        test_matrix(ExactTracking.exact_drift!, M_drift_bmad, L, 10e6, 510998.95; tol=1e-7)

        # GTPSA drift map comparison
        include("bmad_maps/drift.jl")
        L = TPS64{d_z}(1)
        v = transpose(@vars(d_z))
        v = ExactTracking.exact_drift!(1, v, zeros(eltype(v),1,2), L, 10e6, 510998.95)

        coeffs_approx_equal(v_z, v, 1e-8)
        
        # Scalar solenoid jacobian
        L = Float64(1)
        ks = Float64(2)
        test_matrix(ExactTracking.exact_solenoid!, M_solenoid_bmad, L, ks, 10e6, 510998.95; tol=1e-7)

        # GTPSA solenoid jacobian
        L  = TPS64{D}(1)
        ks = TPS64{D}(2)
        test_matrix(ExactTracking.exact_solenoid!, M_solenoid_bmad, L, ks, 10e6, 510998.95; tol=1e-7)

        # GTPSA solenoid map comparison
        include("bmad_maps/solenoid.jl")
        L  = TPS64{d_z}(1)
        ks = TPS64{d_z}(2)
        v = transpose(@vars(d_z))
        v = ExactTracking.exact_solenoid!(1, v, zeros(eltype(v),1,8), L, ks, 10e6, 510998.95)

        coeffs_approx_equal(v_z, v, 1e-8)
    end 
    
end