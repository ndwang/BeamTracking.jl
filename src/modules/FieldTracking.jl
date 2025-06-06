"""
    FieldTracking

Module implementing particle tracking through arbitrary electromagnetic fields using DifferentialEquations.jl.
"""

# Define the Field tracking method
struct Field end

# Number of temporaries needed for a single particle
MAX_TEMPS(::Field) = 0

module FieldTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel
using SciMLBase, OrdinaryDiffEq
const TRACKING_METHOD = Field

# EVOLVE-BLOCK-START
"""
    field_system!(du, u, p, t)

Define the ODE system for particle motion in an electromagnetic field.

# Arguments
- `du`: Vector of derivatives
- `u`: State vector [x, px, y, py, z, pz]
- `p`: Parameters tuple containing (field_func, params)
- `t`: Time variable
"""
function field_system!(du, u, p, t)
    x, px, y, py, z, pz = u
    field_func, params = p
    field = field_func(x, y, z, params)
    
    # Equations of motion
    du[1] = px  # dx/dt = px
    du[2] = field[1]  # dpx/dt = Ex
    du[3] = py  # dy/dt = py
    du[4] = field[2]  # dpy/dt = Ey
    du[5] = pz  # dz/dt = pz
    du[6] = field[3]  # dpz/dt = Ez
end

"""
    field_track!(i, v, work, L, field_func, params, solver)

Track a particle through a drift space with arbitrary field using DifferentialEquations.jl.

# Arguments
- `i`: Particle index
- `v`: Coordinate matrix
- `work`: Work matrix
- `L`: Drift length
- `field_func`: Function that returns the field at a given position (x, y, z)
- `params`: Additional parameters for the field function
- `solver`: ODE solver to use
"""
@makekernel function field_track!(i, v, work, L, field_func, params, solver)
    @inbounds begin
        # Initial state vector [x, px, y, py, z, pz]
        u0 = [v[i,XI], v[i,PXI], v[i,YI], v[i,PYI], v[i,ZI], v[i,PZI]]
        
        # Set up and solve the ODE
        prob = ODEProblem(field_system!, u0, (0.0, L), (field_func, params))
        sol = solve(prob, solver, reltol=1e-8, abstol=1e-8)
        
        # Update final coordinates
        final_state = sol.u[end]
        @FastGTPSA! begin
            v[i,XI] = final_state[1]
            v[i,PXI] = final_state[2]
            v[i,YI] = final_state[3]
            v[i,PYI] = final_state[4]
            v[i,ZI] = final_state[5]
            v[i,PZI] = final_state[6]
        end
    end
    return v
end
# EVOLVE-BLOCK-END

end