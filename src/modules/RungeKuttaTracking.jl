"""
    RungeKuttaTracking

Module implementing particle tracking through arbitrary electromagnetic fields using a 4th order Runge-Kutta method.
"""
module RungeKuttaTracking
using ..BeamTracking, ..StaticArrays
using ..BeamTracking: @makekernel, Coords
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using ..BeamTracking: C_LIGHT, E_CHARGE, vifelse


"""
    multipole_em_field(x, y, z, s, mm, kn, ks)

Compute EM field from multipole moments for RK4 tracking.
Handles ALL multipole orders:
- m=0: solenoid (longitudinal Bz)
- m=1: dipole (transverse By, Bx)
- m≥2: higher-order multipoles (quadrupole, sextupole, etc.)

Returns (Ex, Ey, Ez, Bx, By, Bz) where:
- Bx, By: transverse field from all orders except m=0 (via normalized_field)
- Bz: longitudinal field from m=0 term if present
- Ex, Ey, Ez: zero (static magnetic elements only)
"""
@inline function multipole_em_field(x, y, z, s, mm, kn, ks)
    # Handle empty multipole arrays (pure drift or pure bend)
    if length(mm) == 0
        # No multipole field, return zeros
        bx = zero(x)
        by = zero(y)
        bz = zero(x)
    else
        # Get transverse field components, excluding m=0 (solenoid)
        # normalized_field from multipole.jl handles all orders including m=1 (dipole)
        bx, by = BeamTracking.normalized_field(mm, kn, ks, x, y, 0)

        # Extract longitudinal field from m=0 term if present
        bz = zero(bx)
        if mm[1] == 0
            bz = kn[1]  # Solenoid strength is the m=0 normal component
        end
    end

    # No electric field from static magnets
    return (zero(bx), zero(by), zero(bz), bx, by, bz)
end

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

    # Lorentz force: F = q*(E + v×B)
    E_force_x = charge * Ex
    E_force_y = charge * Ey
    E_force_z = charge * Ez
    B_force_x = charge * (vy*Bz - vz*By)
    B_force_y = charge * (vz*Bx - vx*Bz)
    B_force_z = charge * (vx*By - vy*Bx)

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

    # Total energy for dbeta_ds calculation
    e_tot = p0c * rel_p / beta
    dbeta_ds = mc2^2 * dp_ds * C_LIGHT / e_tot^3

    # Position derivatives: dr/ds = v * dt/ds
    dx_ds = vx * dt_ds
    dy_ds = vy * dt_ds

    # Momentum derivatives: dp_i/ds = F_i * dt/ds / p0c + corrections
    p0 = p0c / C_LIGHT
    dpx_ds = (E_force_x + B_force_x) * dt_ds / p0 + g_bend * pz_p0
    dpy_ds = (E_force_y + B_force_y) * dt_ds / p0

    # Longitudinal coordinate z derivative
    sqrt_1mvt2 = sqrt(1 - vt2_safe)
    dz_ds = rel_dir * (beta / beta_0 - 1) + rel_dir * (sqrt_1mvt2 - 1 - dh_bend) / sqrt_1mvt2 + dbeta_ds * z / beta

    # Energy deviation derivative
    dpz_ds = dp_ds / p0

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
    rk4_step!(v, i, s, h, mm, kn, ks, tracking_params)

Perform a single RK4 step for particle i, updating coordinates in-place.
Uses stack-allocated SVectors for all intermediate values.

# Arguments
- `v`: Coordinate matrix (N_particles × 6)
- `i`: Particle index
- `s`: Current arc length
- `h`: Step size
- `mm`: Multipole orders (StaticArray)
- `kn`: Normal multipole strengths (StaticArray)
- `ks`: Skew multipole strengths (StaticArray)
- `charge`: Particle charge in units of e
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gamsqr_0`: Squared reference Lorentz factor γ₀²
- `g_bend`: Curvature (0 for drift, 1/ρ for bends)
- `p0c`: Reference momentum × c (eV)
- `mc2`: Rest mass energy (eV)
"""
@inline function rk4_step!(v, i, s, h, mm, kn, ks, charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)
    # Extract current state (scalars)
    x = v[i, XI]
    px = v[i, PXI]
    y = v[i, YI]
    py = v[i, PYI]
    z = v[i, ZI]
    pz = v[i, PZI]

    # k1 = f(u, s)
    Ex, Ey, Ez, Bx, By, Bz = multipole_em_field(x, y, z, s, mm, kn, ks)
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
    Ex, Ey, Ez, Bx, By, Bz = multipole_em_field(x2, y2, z2, s + h2, mm, kn, ks)
    k2 = kick_vector(x2, px2, y2, py2, z2, pz2, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # k3 = f(u + h/2 * k2, s + h/2)
    x3 = x + h2 * k2[1]
    px3 = px + h2 * k2[2]
    y3 = y + h2 * k2[3]
    py3 = py + h2 * k2[4]
    z3 = z + h2 * k2[5]
    pz3 = pz + h2 * k2[6]
    Ex, Ey, Ez, Bx, By, Bz = multipole_em_field(x3, y3, z3, s + h2, mm, kn, ks)
    k3 = kick_vector(x3, px3, y3, py3, z3, pz3, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)

    # k4 = f(u + h * k3, s + h)
    x4 = x + h * k3[1]
    px4 = px + h * k3[2]
    y4 = y + h * k3[3]
    py4 = py + h * k3[4]
    z4 = z + h * k3[5]
    pz4 = pz + h * k3[6]
    Ex, Ey, Ez, Bx, By, Bz = multipole_em_field(x4, y4, z4, s + h, mm, kn, ks)
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
    rk4_kernel!(i, coords, beta_0, gamsqr_0, tilde_m, charge, p0c, mc2,
                      s_span, ds_step, g_bend, mm, kn, ks)

Kernelized RK4 tracking through multipole fields.
Compatible with @makekernel and the package's kernel architecture.

The electromagnetic field is computed from multipole moments (mm, kn, ks) using
the multipole_em_field function.
"""
@makekernel function rk4_kernel!(i, coords::Coords, beta_0, gamsqr_0, tilde_m,
                                        charge, p0c, mc2, s_span, ds_step, g_bend, mm, kn, ks)
    # Check if particle is alive at start
    alive_at_start = (coords.state[i] == STATE_ALIVE)

    # Integration loop - only if particle was alive at start
    if alive_at_start
        s_start = s_span[1]
        s_end = s_span[2]
        s = s_start

        v = coords.v
        while s < s_end
            # Calculate remaining distance
            remaining = s_end - s
            # Use ds_step, but if remaining is smaller, use remaining
            h = min(ds_step, remaining)
            # Check momenta before step
            rel_p = 1 + v[i, PZI]
            vt2 = (v[i, PXI] / rel_p)^2 + (v[i, PYI] / rel_p)^2

            # Mark particle as lost if momenta are unphysical (branchless)
            alive = (coords.state[i] == STATE_ALIVE)
            coords.state[i] = vifelse(vt2 >= 1 && alive, STATE_LOST_PZ, coords.state[i])

            # Perform RK4 step
            rk4_step!(v, i, s, h, mm, kn, ks, charge, tilde_m, beta_0, gamsqr_0, g_bend, p0c, mc2)
            s += h
        end
    end
end

end
