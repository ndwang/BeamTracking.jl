"""
  RungeKuttaTracking

Module implementing particle tracking through static magnetic multipole fields
using a fourth-order Runge-Kutta method.
"""
module RungeKuttaTracking
using ..BeamTracking, ..StaticArrays
using ..BeamTracking: @makekernel, Coords
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using ..BeamTracking: C_LIGHT, E_CHARGE, vifelse, normalized_field
using ..BeamTracking: Adapt, FieldMap, fieldmap_em_field


# =====================================================================
# Field Source Structs
# =====================================================================

"""
    MultipoleSource{M,KN,KS,R,GX,GY}

Bundles multipole field-evaluation parameters and reference curvature for the
unified RK4 stepper.
"""
struct MultipoleSource{M,KN,KS,R,GX,GY}
  mm::M
  kn::KN
  ks::KS
  p_over_q_ref::R
  gx::GX
  gy::GY
end

# Keep the original fieldmap-branch constructor for horizontal curvature.
MultipoleSource(mm, kn, ks, p_over_q_ref, gx) =
  MultipoleSource(mm, kn, ks, p_over_q_ref, gx, zero(gx))

"""
    FieldMapSource{FM,T,GX,GY}

Bundles a field map, its longitudinal offset, and reference curvature for the
unified RK4 stepper.
"""
struct FieldMapSource{FM,T,GX,GY}
  fieldmap::FM
  z_offset::T
  gx::GX
  gy::GY
end

# Field maps are normally attached to straight elements. This constructor
# preserves the original three-argument API while allowing tilted curvature in
# the full constructor.
FieldMapSource(fieldmap, z_offset, gx) =
  FieldMapSource(fieldmap, z_offset, gx, zero(gx))

Adapt.@adapt_structure MultipoleSource
Adapt.@adapt_structure FieldMapSource

# Kernel arguments are lowered before launch and evaluated for each batch or
# particle time. Preserve that behavior when parameters are wrapped in a field
# source instead of passed as separate tuple entries.
BeamTracking.batch_lower(f::MultipoleSource) = MultipoleSource(
  BeamTracking.batch_lower(f.mm),
  BeamTracking.batch_lower(f.kn),
  BeamTracking.batch_lower(f.ks),
  BeamTracking.batch_lower(f.p_over_q_ref),
  BeamTracking.batch_lower(f.gx),
  BeamTracking.batch_lower(f.gy),
)
BeamTracking.static_batchcheck(f::MultipoleSource) =
  BeamTracking.static_batchcheck((f.mm, f.kn, f.ks, f.p_over_q_ref, f.gx, f.gy))
BeamTracking.beval(f::MultipoleSource, i) = MultipoleSource(
  BeamTracking.beval(f.mm, i),
  BeamTracking.beval(f.kn, i),
  BeamTracking.beval(f.ks, i),
  BeamTracking.beval(f.p_over_q_ref, i),
  BeamTracking.beval(f.gx, i),
  BeamTracking.beval(f.gy, i),
)

BeamTracking.time_lower(f::MultipoleSource) = MultipoleSource(
  BeamTracking.time_lower(f.mm),
  BeamTracking.time_lower(f.kn),
  BeamTracking.time_lower(f.ks),
  BeamTracking.time_lower(f.p_over_q_ref),
  BeamTracking.time_lower(f.gx),
  BeamTracking.time_lower(f.gy),
)
BeamTracking.static_timecheck(f::MultipoleSource) =
  BeamTracking.static_timecheck((f.mm, f.kn, f.ks, f.p_over_q_ref, f.gx, f.gy))
BeamTracking.teval(f::MultipoleSource, t) = MultipoleSource(
  BeamTracking.teval(f.mm, t),
  BeamTracking.teval(f.kn, t),
  BeamTracking.teval(f.ks, t),
  BeamTracking.teval(f.p_over_q_ref, t),
  BeamTracking.teval(f.gx, t),
  BeamTracking.teval(f.gy, t),
)

BeamTracking.batch_lower(f::FieldMapSource) = FieldMapSource(
  BeamTracking.batch_lower(f.fieldmap),
  BeamTracking.batch_lower(f.z_offset),
  BeamTracking.batch_lower(f.gx),
  BeamTracking.batch_lower(f.gy),
)
BeamTracking.static_batchcheck(f::FieldMapSource) =
  BeamTracking.static_batchcheck((f.fieldmap, f.z_offset, f.gx, f.gy))
BeamTracking.beval(f::FieldMapSource, i) = FieldMapSource(
  BeamTracking.beval(f.fieldmap, i),
  BeamTracking.beval(f.z_offset, i),
  BeamTracking.beval(f.gx, i),
  BeamTracking.beval(f.gy, i),
)

BeamTracking.time_lower(f::FieldMapSource) = FieldMapSource(
  BeamTracking.time_lower(f.fieldmap),
  BeamTracking.time_lower(f.z_offset),
  BeamTracking.time_lower(f.gx),
  BeamTracking.time_lower(f.gy),
)
BeamTracking.static_timecheck(f::FieldMapSource) =
  BeamTracking.static_timecheck((f.fieldmap, f.z_offset, f.gx, f.gy))
BeamTracking.teval(f::FieldMapSource, t) = FieldMapSource(
  BeamTracking.teval(f.fieldmap, t),
  BeamTracking.teval(f.z_offset, t),
  BeamTracking.teval(f.gx, t),
  BeamTracking.teval(f.gy, t),
)

@inline eval_em_field(f::MultipoleSource, x, y, z, pz, s) =
  multipole_em_field(x, y, z, s, f.mm, f.kn, f.ks, f.p_over_q_ref)

@inline eval_em_field(f::FieldMapSource, x, y, z, pz, s) =
  fieldmap_em_field(x, y, z, pz, s, f.fieldmap, f.z_offset)

@inline get_curvature(f) = (f.gx, f.gy)
@inline get_g_bend(f) = f.gx


"""
  multipole_em_field(x, y, z, s, mm, kn, ks, p_over_q_ref)

Compute EM field from multipole moments for RK4 tracking.
Handles ALL multipole orders:
- m=0: solenoid (longitudinal Bz)
- m=1: dipole (transverse By, Bx)
- m≥2: higher-order multipoles (quadrupole, sextupole, etc.)

Returns (Ex, Ey, Ez, Bx, By, Bz) in physical units (Tesla for B, V/m for E) where:
- Bx, By: transverse field from all orders except m=0
- Bz: longitudinal field from m=0 term if present
- Ex, Ey, Ez: zero (static magnetic elements only)
"""
@inline function multipole_em_field(x, y, z, s, mm::SVector{0}, kn, ks, p_over_q_ref)
  return (zero(x), zero(x), zero(x), zero(x), zero(x), zero(x))
end

@inline function multipole_em_field(x, y, z, s, mm::SVector{N}, kn, ks, p_over_q_ref) where N
  bx, by = normalized_field(mm, kn, ks, x, y, 0)
  is_solenoid = (mm[1] == 0)
  bz = vifelse(is_solenoid, kn[1], zero(x))

  # Convert from normalized (field/Bρ) to physical units (Tesla)
  return (zero(x), zero(x), zero(x), bx * p_over_q_ref, by * p_over_q_ref, bz * p_over_q_ref)
end

"""
  kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
        charge, tilde_m, beta_0, gx, gy, p0c, mc2)

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
- `gx`, `gy`: Horizontal and vertical reference curvature components
- `p0c`: Reference momentum × c (eV)
- `mc2`: Rest mass energy (eV)
"""
@inline function kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)
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
  inv_gamma_v = sqrt(rel_p2 + tilde_m^2)
  beta = rel_p / inv_gamma_v
  
  inv_beta_c = 1 / (beta * C_LIGHT)

  # Longitudinal velocity component
  rel_dir = 1  # +1 for forward tracking
  vz_norm = sqrt(1 - vt2_safe) * rel_dir
  vx = beta * C_LIGHT * vt_x
  vy = beta * C_LIGHT * vt_y
  vz = beta * C_LIGHT * vz_norm

  # Combine charge and reference momentum once, then apply the Lorentz force.
  qp0 = charge * C_LIGHT / p0c
  Fx = Ex + (vy*Bz - vz*By)
  Fy = Ey + (vz*Bx - vx*Bz)

  # Time derivative w.r.t. arc length
  dh_bend = x * gx + y * gy  # Longitudinal distance deviation
  abs_vz = abs(vz)
  abs_vz_safe = vifelse(good_momenta, abs_vz, one(abs_vz))  # Avoid division by zero
  dt_ds = rel_dir * (1 + dh_bend) / abs_vz_safe

  # Longitudinal momentum (normalized)
  pz_p0 = rel_p * rel_dir * abs_vz * inv_beta_c

  # Magnetic fields do no work, so only the electric field changes energy.
  E_dot_v = Ex*vx + Ey*vy + Ez*vz
  dpz_ds = qp0 * E_dot_v * dt_ds * inv_beta_c

  # Total energy for dbeta_ds calculation
  e_tot = p0c * rel_p / beta
  dbeta_ds = mc2^2 * dpz_ds * p0c / e_tot^3

  # Position derivatives: dr/ds = v * dt/ds
  dx_ds = vx * dt_ds
  dy_ds = vy * dt_ds

  # Momentum derivatives, including both components of reference curvature.
  dpx_ds = qp0 * Fx * dt_ds + gx * pz_p0
  dpy_ds = qp0 * Fy * dt_ds + gy * pz_p0

  # Longitudinal coordinate z derivative
  sqrt_1mvt2 = sqrt(1 - vt2_safe)
  dz_ds = rel_dir * (beta / beta_0 - 1) + rel_dir * (sqrt_1mvt2 - 1 - dh_bend) / sqrt_1mvt2 + dbeta_ds * z / beta

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
  rk4_step!(coords, i, s, h, field, charge, tilde_m, beta_0, p0c, mc2)

Perform a single RK4 step for particle i, updating coordinates in-place.
Only updates state if particle is alive.

The `field` argument is a `MultipoleSource` or `FieldMapSource` and carries
both the field-evaluation data and reference curvature.
"""
@inline function rk4_step!(coords, i, s, h, field, charge, tilde_m, beta_0, p0c, mc2)
  gx, gy = get_curvature(field)

  # Check if particle is alive
  alive = (coords.state[i] == STATE_ALIVE)
  
  # Extract current particle
  v = coords.v
  x = v[i, XI]
  px = v[i, PXI]
  y = v[i, YI]
  py = v[i, PYI]
  z = v[i, ZI]
  pz = v[i, PZI]

  # k1 = f(u, s)
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field(field, x, y, z, pz, s)
  k1 = kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

  # k2 = f(u + h/2 * k1, s + h/2)
  h2 = h / 2
  x2 = x + h2 * k1[1]
  px2 = px + h2 * k1[2]
  y2 = y + h2 * k1[3]
  py2 = py + h2 * k1[4]
  z2 = z + h2 * k1[5]
  pz2 = pz + h2 * k1[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field(field, x2, y2, z2, pz2, s + h2)
  k2 = kick_vector(x2, px2, y2, py2, z2, pz2, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

  # k3 = f(u + h/2 * k2, s + h/2)
  x3 = x + h2 * k2[1]
  px3 = px + h2 * k2[2]
  y3 = y + h2 * k2[3]
  py3 = py + h2 * k2[4]
  z3 = z + h2 * k2[5]
  pz3 = pz + h2 * k2[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field(field, x3, y3, z3, pz3, s + h2)
  k3 = kick_vector(x3, px3, y3, py3, z3, pz3, s + h2, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

  # k4 = f(u + h * k3, s + h)
  x4 = x + h * k3[1]
  px4 = px + h * k3[2]
  y4 = y + h * k3[3]
  py4 = py + h * k3[4]
  z4 = z + h * k3[5]
  pz4 = pz + h * k3[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field(field, x4, y4, z4, pz4, s + h)
  k4 = kick_vector(x4, px4, y4, py4, z4, pz4, s + h, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

  # Update state: u += h/6 * (k1 + 2*k2 + 2*k3 + k4)
  # Only update if particle is alive
  h6 = h / 6
  v[i, XI] = vifelse(alive, x + h6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]), v[i, XI])
  v[i, PXI] = vifelse(alive, px + h6 * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]), v[i, PXI])
  v[i, YI] = vifelse(alive, y + h6 * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]), v[i, YI])
  v[i, PYI] = vifelse(alive, py + h6 * (k1[4] + 2*k2[4] + 2*k3[4] + k4[4]), v[i, PYI])
  v[i, ZI] = vifelse(alive, z + h6 * (k1[5] + 2*k2[5] + 2*k3[5] + k4[5]), v[i, ZI])
  v[i, PZI] = vifelse(alive, pz + h6 * (k1[6] + 2*k2[6] + 2*k3[6] + k4[6]), v[i, PZI])
end

"""
  rk4_kernel!(i, coords, beta_0, tilde_m, charge, p0c, mc2,
              L, ds_step, n_steps, field)

Kernelized RK4 tracking through electromagnetic fields.
Compatible with @makekernel and the package's kernel architecture.

The `field` argument selects analytical multipole or gridded field-map
evaluation through dispatch.
"""
@makekernel function rk4_kernel!(i, coords::Coords, beta_0, tilde_m,
                                charge, p0c, mc2, L, ds_step, n_steps, field)
  v = coords.v
  
  for step in 1:n_steps
    # Derive the position from the step index instead of accumulating it.
    # Accumulation can place the final RK evaluation just beyond a field-map
    # boundary because of floating-point roundoff.
    s = (step - 1) * ds_step

    # Check if particle is lost
    rel_p = 1 + v[i, PZI]
    inv_rel_p = 1 / rel_p
    vt2 = (v[i, PXI] * inv_rel_p)^2 + (v[i, PYI] * inv_rel_p)^2
    alive = (coords.state[i] == STATE_ALIVE)
    # Mark particle as lost
    coords.state[i] = vifelse((vt2 >= 1) & alive, STATE_LOST_PZ, coords.state[i])

    # Perform RK4 step (check for alive status is now inside rk4_step!)
    rk4_step!(coords, i, s, ds_step, field, charge, tilde_m, beta_0, p0c, mc2)

    # The common path performs the final callback after exit processing.
    if step != n_steps
      next_s = step * ds_step
      BeamTracking.execute_callbacks(i, coords, next_s, next_s / (beta_0 * C_LIGHT))
    end
  end
end

end
