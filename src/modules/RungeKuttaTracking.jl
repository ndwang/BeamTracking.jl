struct RungeKutta end

"""
    RungeKuttaFieldTracking

Module implementing particle tracking through arbitrary electromagnetic fields using a 4th order Runge-Kutta method.
"""
module RungeKuttaTracking
using ..BeamTracking, ..StaticArrays
using ..BeamTracking: @makekernel, Coords
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST, STATE_LOST_PZ
using ..BeamTracking: C_LIGHT, E_CHARGE, vifelse

const TRACKING_METHOD = RungeKutta

"""
    kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

Calculate the derivative vector du/ds for relativistic particle tracking.
Returns an SVector{6} containing [dx/ds, dpx/ds, dy/ds, dpy/ds, dz/ds, dpz/ds].

Uses branchless operations for GPU/SIMD compatibility. For unphysical momenta,
returns zero derivatives (caller should mark particle as lost).

# Arguments
- `x, px, y, py, z, pz`: State vector components
- `s`: Arc length position
- `Ex, Ey, Ez`: Electric field components (V/m)
- `Bx, By, Bz`: Magnetic field components (T)
- `charge`: Particle charge in units of e
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gamsqr_0`: Squared reference Lorentz factor γ₀²
- `g_bend`: Curvature (0 for drift, 1/ρ for bends)
- `p0c`: Reference momentum × c (eV)
- `mc2`: Rest mass energy (eV)
"""
@inline function kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                              charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)
    # Relative momentum
    rel_p = 1 + pz

    # Transverse velocity components (normalized)
    vt_x = px / rel_p
    vt_y = py / rel_p
    vt2 = vt_x^2 + vt_y^2

    # Check for unphysical momenta (branchless)
    vt2_1 = one(vt2)
    good_momenta = (vt2 < vt2_1)
    vt2_safe = vifelse(good_momenta, vt2, zero(vt2))

    # Particle beta and velocity
    rel_p2 = rel_p^2
    beta = rel_p / sqrt(rel_p2 + tilde_m^2)

    # Longitudinal velocity component
    rel_dir = 1  # +1 for forward tracking
    vz_norm = sqrt(1 - vt2_safe) * rel_dir
    vx = beta * C_LIGHT * vt_x
    vy = beta * C_LIGHT * vt_y
    vz = beta * C_LIGHT * vz_norm

    # Lorentz force: F = q*(E + v×B) in Newtons
    # E_force in eV/m = (charge * e) * E[V/m]
    E_force_x = charge * E_CHARGE * Ex
    E_force_y = charge * E_CHARGE * Ey
    E_force_z = charge * E_CHARGE * Ez

    # B_force = q*(v × B), component by component
    B_force_x = charge * E_CHARGE * (vy*Bz - vz*By)
    B_force_y = charge * E_CHARGE * (vz*Bx - vx*Bz)
    B_force_z = charge * E_CHARGE * (vx*By - vy*Bx)

    # Time derivative w.r.t. arc length
    dh_bend = x * g_bend  # Longitudinal distance deviation
    abs_vz = abs(vz)
    abs_vz_safe = vifelse(good_momenta, abs_vz, one(abs_vz))  # Avoid division by zero
    dt_ds = rel_dir * (1 + dh_bend) / abs_vz_safe

    # Longitudinal momentum (normalized)
    pz_p0 = rel_p * rel_dir * abs_vz / (beta * C_LIGHT)

    # Energy derivative: dp/ds = (F · v) * dt/ds / (β*c)
    F_dot_v = E_force_x*vx + E_force_y*vy + E_force_z*vz
    dp_ds = F_dot_v * dt_ds / (beta * C_LIGHT)

    # Position derivatives: dr/ds = v * dt/ds
    dx_ds = vx * dt_ds
    dy_ds = vy * dt_ds

    # Momentum derivatives: dp_i/ds = F_i * dt/ds / p0c + corrections
    dpx_ds = (E_force_x + B_force_x) * dt_ds / p0c + g_bend * pz_p0
    dpy_ds = (E_force_y + B_force_y) * dt_ds / p0c

    # Longitudinal coordinate z derivative
    # Simplified formula (can add full Fortran version later)
    sqrt_1mvt2 = sqrt(1 - vt2_safe)
    dz_ds = rel_dir * (beta / beta_0 - 1) + rel_dir * (sqrt_1mvt2 - dh_bend) / sqrt_1mvt2

    # Energy deviation derivative
    dpz_ds = dp_ds / p0c

    # Return zero derivatives if momenta are unphysical (branchless)
    zero_deriv = zero(dx_ds)
    return SVector(
        vifelse(good_momenta, dx_ds, zero_deriv),
        vifelse(good_momenta, dpx_ds, zero_deriv),
        vifelse(good_momenta, dy_ds, zero_deriv),
        vifelse(good_momenta, dpy_ds, zero_deriv),
        vifelse(good_momenta, dz_ds, zero_deriv),
        vifelse(good_momenta, dpz_ds, zero_deriv)
    )
end

"""
    rk4_step!(v, i, s, h, field_func, field_params, tracking_params)

Perform a single RK4 step for particle i, updating coordinates in-place.
Uses stack-allocated SVectors for all intermediate values.

# Arguments
- `v`: Coordinate matrix (N_particles × 6)
- `i`: Particle index
- `s`: Current arc length
- `h`: Step size
- `field_func`: Function returning (Ex, Ey, Ez, Bx, By, Bz) = field_func(x, px, y, py, z, pz, s, field_params)
- `field_params`: Parameters for field function
- `tracking_params`: Tuple of (charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)
"""
@inline function rk4_step!(v, i, s, h, field_func, field_params, tracking_params)
    # Unpack tracking parameters
    charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2 = tracking_params

    # Extract current state (scalars)
    x = v[i, XI]
    px = v[i, PXI]
    y = v[i, YI]
    py = v[i, PYI]
    z = v[i, ZI]
    pz = v[i, PZI]

    # k1 = f(u, s)
    Ex, Ey, Ez, Bx, By, Bz = field_func(x, px, y, py, z, pz, s, field_params)
    k1 = kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # k2 = f(u + h/2 * k1, s + h/2)
    h2 = h / 2
    x2 = x + h2 * k1[1]
    px2 = px + h2 * k1[2]
    y2 = y + h2 * k1[3]
    py2 = py + h2 * k1[4]
    z2 = z + h2 * k1[5]
    pz2 = pz + h2 * k1[6]
    Ex, Ey, Ez, Bx, By, Bz = field_func(x2, px2, y2, py2, z2, pz2, s + h2, field_params)
    k2 = kick_vector(x2, px2, y2, py2, z2, pz2, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # k3 = f(u + h/2 * k2, s + h/2)
    x3 = x + h2 * k2[1]
    px3 = px + h2 * k2[2]
    y3 = y + h2 * k2[3]
    py3 = py + h2 * k2[4]
    z3 = z + h2 * k2[5]
    pz3 = pz + h2 * k2[6]
    Ex, Ey, Ez, Bx, By, Bz = field_func(x3, px3, y3, py3, z3, pz3, s + h2, field_params)
    k3 = kick_vector(x3, px3, y3, py3, z3, pz3, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # k4 = f(u + h * k3, s + h)
    x4 = x + h * k3[1]
    px4 = px + h * k3[2]
    y4 = y + h * k3[3]
    py4 = py + h * k3[4]
    z4 = z + h * k3[5]
    pz4 = pz + h * k3[6]
    Ex, Ey, Ez, Bx, By, Bz = field_func(x4, px4, y4, py4, z4, pz4, s + h, field_params)
    k4 = kick_vector(x4, px4, y4, py4, z4, pz4, s + h, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # Update state: u += h/6 * (k1 + 2*k2 + 2*k3 + k4)
    h6 = h / 6
    v[i, XI] = x + h6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    v[i, PXI] = px + h6 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    v[i, YI] = y + h6 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])
    v[i, PYI] = py + h6 * (k1[4] + 2*k2[4] + 2*k3[4] + k4[4])
    v[i, ZI] = z + h6 * (k1[5] + 2*k2[5] + 2*k3[5] + k4[5])
    v[i, PZI] = pz + h6 * (k1[6] + 2*k2[6] + 2*k3[6] + k4[6])
end

"""
    rk4_track!(i, coords, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2,
               s_span, field_func, field_params, n_steps, g_bend)

Track particle through arbitrary electromagnetic fields using RK4.

# Arguments
- `i`: Particle index
- `coords`: Coords object containing particle coordinates and state
- `beta_0`: Reference velocity β₀ = v₀/c
- `gamsqr_0`: Squared reference Lorentz factor γ₀²
- `tilde_m`: Normalized mass mc²/(p₀c)
- `charge`: Particle charge in units of e
- `p0c`: Reference momentum × c (eV)
- `mc2`: Rest mass energy (eV)
- `s_span`: Arc length span [s_start, s_end]
- `field_func`: Function returning (Ex, Ey, Ez, Bx, By, Bz) = field_func(x, px, y, py, z, pz, s, field_params)
- `field_params`: Parameters for field function
- `n_steps`: Number of integration steps
- `g_bend`: Curvature (0 for drift, 1/ρ for bends)
"""
function rk4_track!(i, coords::Coords, beta_0, gamsqr_0, tilde_m,
                    charge, p0c, mc2, s_span, field_func,
                    field_params, n_steps, g_bend)
    v = coords.v

    # Check if particle is alive at start
    alive_at_start = (coords.state[i] == STATE_ALIVE)

    # Skip if particle is not alive
    if !alive_at_start
        return
    end

    tracking_params = (charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # Integration
    h = (s_span[2] - s_span[1]) / n_steps
    s = s_span[1]

    for step in 1:n_steps
        # Check momenta before step
        rel_p = 1 + v[i, PZI]
        vt2 = (v[i, PXI] / rel_p)^2 + (v[i, PYI] / rel_p)^2
        good_momenta = (vt2 < 1)

        # Mark particle as lost if momenta are unphysical
        alive = (coords.state[i] == STATE_ALIVE)
        coords.state[i] = vifelse(!good_momenta & alive, STATE_LOST_PZ, coords.state[i])

        # Perform RK4 step (will return zero derivatives if lost)
        rk4_step!(v, i, s, h, field_func, field_params, tracking_params)
        s += h
    end
end

end