using Test
using BeamTracking
using StaticArrays

# Test field_system! with a uniform electric field
@testset "FieldTracking" begin
    # Define a simple uniform electric field in x-direction
    function uniform_field(x, y, z, params)
        return SVector(1.0, 0.0, 0.0)
    end

    # Test initial conditions
    du = zeros(6)
    u = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Initial position at x=1, rest at origin
    p = (uniform_field, nothing)
    t = 0.0

    # Call field_system!
    FieldTracking.field_system!(du, u, p, t)

    # Test derivatives
    @test du[1] ≈ 0.0  # dx/dt = px = 0
    @test du[2] ≈ 1.0  # dpx/dt = Ex = 1
    @test du[3] ≈ 0.0  # dy/dt = py = 0
    @test du[4] ≈ 0.0  # dpy/dt = Ey = 0
    @test du[5] ≈ 0.0  # dz/dt = pz = 0
    @test du[6] ≈ 0.0  # dpz/dt = Ez = 0

    # Test with non-zero initial momentum
    u = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]  # Initial momentum px=1
    du = zeros(6)
    BeamTracking.field_system!(du, u, p, t)

    @test du[1] ≈ 1.0  # dx/dt = px = 1
    @test du[2] ≈ 1.0  # dpx/dt = Ex = 1
    @test du[3] ≈ 0.0  # dy/dt = py = 0
    @test du[4] ≈ 0.0  # dpy/dt = Ey = 0
    @test du[5] ≈ 0.0  # dz/dt = pz = 0
    @test du[6] ≈ 0.0  # dpz/dt = Ez = 0
end 