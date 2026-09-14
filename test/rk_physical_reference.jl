# Physical-unit implementation from cf7e2d1, retained for numerical regression.
module RKPhysicalReference
using BeamTracking, StaticArrays
using BeamTracking: C_LIGHT, EMField, vifelse

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

  # Lorentz force: F = q*(E + v×B)
  E_force_x = charge * Ex
  E_force_y = charge * Ey
  E_force_z = charge * Ez
  B_force_x = charge * (vy*Bz - vz*By)
  B_force_y = charge * (vz*Bx - vx*Bz)
  B_force_z = charge * (vx*By - vy*Bx)

  # Time derivative w.r.t. arc length
  dh_bend = x * gx + y * gy  # Longitudinal distance deviation
  abs_vz = abs(vz)
  abs_vz_safe = vifelse(good_momenta, abs_vz, one(abs_vz))  # Avoid division by zero
  dt_ds = rel_dir * (1 + dh_bend) / abs_vz_safe

  # Longitudinal momentum (normalized)
  pz_p0 = rel_p * rel_dir * abs_vz * inv_beta_c

  # Energy derivative: dp/ds = (F · v) * dt/ds * inv_beta_c
  F_dot_v = E_force_x*vx + E_force_y*vy + E_force_z*vz
  dp_ds = F_dot_v * dt_ds * inv_beta_c

  # Total energy for dbeta_ds calculation
  e_tot = p0c * rel_p / beta
  dbeta_ds = mc2^2 * dp_ds * C_LIGHT / e_tot^3

  # Position derivatives: dr/ds = v * dt/ds
  dx_ds = vx * dt_ds
  dy_ds = vy * dt_ds

  # Momentum derivatives: dp_i/ds = F_i * dt/ds / p0c + corrections
  p0 = p0c / C_LIGHT
  dpx_ds = (E_force_x + B_force_x) * dt_ds / p0 + gx * pz_p0
  dpy_ds = (E_force_y + B_force_y) * dt_ds / p0 + gy * pz_p0

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

@inline function kick_vector(x, px, y, py, z, pz, s, field::EMField,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  Ex, Ey, Ez = field.E
  Bx, By, Bz = field.B
  return kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                     charge, tilde_m, beta_0, gx, gy, p0c, mc2)
end

end
