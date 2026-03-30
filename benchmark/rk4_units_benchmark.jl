"""
Benchmark: Physical units vs normalized units in RK4 kernel.

Tests a 2×2 matrix:
  - Field source: analytical multipole vs field map (3D rectangular grid)
  - Unit convention: physical (current) vs normalized (proposed)

On both CPU and GPU (CUDA).

Uses the standard BeamTracking kernel launch sequence (KernelChain + launch!).
"""

using BeamTracking
using BeamTracking: Species, massof, chargeof, R_to_beta_gamma, R_to_pc, pc_to_R,
                    RungeKuttaTracking, Bunch, Coords, STATE_ALIVE,
                    KernelCall, KernelChain, launch!, @makekernel,
                    XI, PXI, YI, PYI, ZI, PZI, STATE_LOST_PZ,
                    C_LIGHT, E_CHARGE, vifelse, normalized_field,
                    FieldMap, RectGrid3D, fieldmap_em_field
using BeamTracking.RungeKuttaTracking: MultipoleSource, FieldMapSource, kick_vector
using StaticArrays
using BenchmarkTools
using Statistics: std
using BeamTracking: Adapt

_has_cuda = false
try
  using CUDA
  global _has_cuda = CUDA.functional()
catch
end

# ============================================================
# Physical-units kernels: already in RungeKuttaTracking.rk4_kernel!
# Works with both MultipoleSource and FieldMapSource via eval_em_field dispatch
# ============================================================

# ============================================================
# Superposition kernels: multipole + field map evaluated and added
# ============================================================
module SuperpositionPhysical
using BeamTracking: @makekernel, Coords
using BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using BeamTracking: C_LIGHT, vifelse
using BeamTracking.RungeKuttaTracking: MultipoleSource, FieldMapSource, eval_em_field, get_g_bend, kick_vector
using StaticArrays

@inline function rk4_step_super!(coords, i, s, h, mp_field, fm_field, charge, tilde_m, beta_0, p0c, mc2)
  g_bend = get_g_bend(mp_field)
  alive = (coords.state[i] == STATE_ALIVE)
  v = coords.v
  x  = v[i, XI]; px = v[i, PXI]
  y  = v[i, YI]; py = v[i, PYI]
  z  = v[i, ZI]; pz = v[i, PZI]

  # Evaluate both sources in physical units, add
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_em_field(mp_field, x, y, z, pz, s)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_em_field(fm_field, x, y, z, pz, s)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k1 = kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, g_bend, p0c, mc2)

  h2 = h / 2
  x2 = x+h2*k1[1]; px2 = px+h2*k1[2]; y2 = y+h2*k1[3]
  py2 = py+h2*k1[4]; z2 = z+h2*k1[5]; pz2 = pz+h2*k1[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_em_field(mp_field, x2, y2, z2, pz2, s+h2)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_em_field(fm_field, x2, y2, z2, pz2, s+h2)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k2 = kick_vector(x2, px2, y2, py2, z2, pz2, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, g_bend, p0c, mc2)

  x3 = x+h2*k2[1]; px3 = px+h2*k2[2]; y3 = y+h2*k2[3]
  py3 = py+h2*k2[4]; z3 = z+h2*k2[5]; pz3 = pz+h2*k2[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_em_field(mp_field, x3, y3, z3, pz3, s+h2)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_em_field(fm_field, x3, y3, z3, pz3, s+h2)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k3 = kick_vector(x3, px3, y3, py3, z3, pz3, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, g_bend, p0c, mc2)

  x4 = x+h*k3[1]; px4 = px+h*k3[2]; y4 = y+h*k3[3]
  py4 = py+h*k3[4]; z4 = z+h*k3[5]; pz4 = pz+h*k3[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_em_field(mp_field, x4, y4, z4, pz4, s+h)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_em_field(fm_field, x4, y4, z4, pz4, s+h)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k4 = kick_vector(x4, px4, y4, py4, z4, pz4, s+h, Ex, Ey, Ez, Bx, By, Bz,
                charge, tilde_m, beta_0, g_bend, p0c, mc2)

  h6 = h / 6
  v[i, XI]  = vifelse(alive, x  + h6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]), v[i, XI])
  v[i, PXI] = vifelse(alive, px + h6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]), v[i, PXI])
  v[i, YI]  = vifelse(alive, y  + h6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]), v[i, YI])
  v[i, PYI] = vifelse(alive, py + h6*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]), v[i, PYI])
  v[i, ZI]  = vifelse(alive, z  + h6*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5]), v[i, ZI])
  v[i, PZI] = vifelse(alive, pz + h6*(k1[6] + 2*k2[6] + 2*k3[6] + k4[6]), v[i, PZI])
end

@makekernel function rk4_kernel_super!(i, coords::Coords, beta_0, tilde_m,
                                charge, p0c, mc2, s_span, ds_step, mp_field, fm_field)
  s_start = s_span[1]; s_end = s_span[2]; s = s_start
  v = coords.v
  total_distance = s_end - s_start
  n_steps = ceil(Int, total_distance / ds_step)
  for step in 1:n_steps
    remaining = s_end - s; h = min(ds_step, remaining)
    rel_p = 1 + v[i, PZI]; inv_rel_p = 1 / rel_p
    vt2 = (v[i, PXI] * inv_rel_p)^2 + (v[i, PYI] * inv_rel_p)^2
    alive = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse((vt2 >= 1) & alive, STATE_LOST_PZ, coords.state[i])
    rk4_step_super!(coords, i, s, h, mp_field, fm_field, charge, tilde_m, beta_0, p0c, mc2)
    s += h
  end
end
end # SuperpositionPhysical

module SuperpositionNormalized
using BeamTracking: @makekernel, Coords
using BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using BeamTracking: C_LIGHT, vifelse, normalized_field, fieldmap_em_field
using BeamTracking.RungeKuttaTracking: MultipoleSource, FieldMapSource
using StaticArrays

# Multipole: return normalized (no * Bρ)
@inline function eval_mp_normalized(f::MultipoleSource, x, y, z, pz, s)
  mm = f.mm; kn = f.kn; ks = f.ks
  bx, by = normalized_field(mm, kn, ks, x, y, 0)
  is_solenoid = (mm[1] == 0)
  bz = vifelse(is_solenoid, kn[1], zero(x))
  return (zero(x), zero(x), zero(x), bx, by, bz)
end

# Field map: return physical, then divide by Bρ to normalize
@inline function eval_fm_normalized(f::FieldMapSource, x, y, z, pz, s, inv_brho)
  Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(x, y, z, pz, s, f.fieldmap, f.z_offset)
  return (Ex*inv_brho, Ey*inv_brho, Ez*inv_brho,
          Bx*inv_brho, By*inv_brho, Bz*inv_brho)
end

# Simplified kick_vector (same as NormalizedKernel)
@inline function kick_vector_normalized(x, px, y, py, z, pz, s, en_x, en_y, en_z, bn_x, bn_y, bn_z,
                tilde_m, beta_0, g_bend, p0c, mc2)
  rel_p = 1 + pz
  vt_x = px / rel_p
  vt_y = py / rel_p
  vt2 = vt_x^2 + vt_y^2
  vt2_1 = one(vt2)
  good_momenta = (vt2 < vt2_1)
  vt2_safe = vifelse(good_momenta, vt2, zero(vt2))
  rel_p2 = rel_p^2
  inv_gamma_v = sqrt(rel_p2 + tilde_m^2)
  beta = rel_p / inv_gamma_v
  inv_beta_c = 1 / (beta * C_LIGHT)
  rel_dir = 1
  vz_norm = sqrt(1 - vt2_safe) * rel_dir
  vx = beta * C_LIGHT * vt_x
  vy = beta * C_LIGHT * vt_y
  vz = beta * C_LIGHT * vz_norm
  Fn_x = en_x + (vy*bn_z - vz*bn_y)
  Fn_y = en_y + (vz*bn_x - vx*bn_z)
  dh_bend = x * g_bend
  abs_vz = abs(vz)
  abs_vz_safe = vifelse(good_momenta, abs_vz, one(abs_vz))
  dt_ds = rel_dir * (1 + dh_bend) / abs_vz_safe
  pz_p0 = rel_p * rel_dir * abs_vz * inv_beta_c
  en_dot_v = en_x*vx + en_y*vy + en_z*vz
  dpz_ds_val = en_dot_v * dt_ds * inv_beta_c
  e_tot = p0c * rel_p / beta
  dbeta_ds = mc2^2 * dpz_ds_val * p0c / e_tot^3
  dx_ds = vx * dt_ds
  dy_ds = vy * dt_ds
  dpx_ds = Fn_x * dt_ds + g_bend * pz_p0
  dpy_ds = Fn_y * dt_ds
  sqrt_1mvt2 = sqrt(1 - vt2_safe)
  dz_ds = rel_dir * (beta / beta_0 - 1) + rel_dir * (sqrt_1mvt2 - 1 - dh_bend) / sqrt_1mvt2 + dbeta_ds * z / beta
  zero_deriv = zero(dx_ds)
  return SVector(
    vifelse(good_momenta, dx_ds, zero_deriv),
    vifelse(good_momenta, dpx_ds, zero_deriv),
    vifelse(good_momenta, dy_ds, zero_deriv),
    vifelse(good_momenta, dpy_ds, zero_deriv),
    vifelse(good_momenta, dz_ds, zero_deriv),
    vifelse(good_momenta, dpz_ds_val, zero_deriv)
  )
end

@inline function rk4_step_super_norm!(coords, i, s, h, mp_field, fm_field, tilde_m, beta_0, p0c, mc2, inv_brho)
  g_bend = mp_field.g_bend
  alive = (coords.state[i] == STATE_ALIVE)
  v = coords.v
  x  = v[i, XI]; px = v[i, PXI]
  y  = v[i, YI]; py = v[i, PYI]
  z  = v[i, ZI]; pz = v[i, PZI]

  # Evaluate both in normalized, add
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_mp_normalized(mp_field, x, y, z, pz, s)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_fm_normalized(fm_field, x, y, z, pz, s, inv_brho)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k1 = kick_vector_normalized(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  h2 = h / 2
  x2 = x+h2*k1[1]; px2 = px+h2*k1[2]; y2 = y+h2*k1[3]
  py2 = py+h2*k1[4]; z2 = z+h2*k1[5]; pz2 = pz+h2*k1[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_mp_normalized(mp_field, x2, y2, z2, pz2, s+h2)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_fm_normalized(fm_field, x2, y2, z2, pz2, s+h2, inv_brho)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k2 = kick_vector_normalized(x2, px2, y2, py2, z2, pz2, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  x3 = x+h2*k2[1]; px3 = px+h2*k2[2]; y3 = y+h2*k2[3]
  py3 = py+h2*k2[4]; z3 = z+h2*k2[5]; pz3 = pz+h2*k2[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_mp_normalized(mp_field, x3, y3, z3, pz3, s+h2)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_fm_normalized(fm_field, x3, y3, z3, pz3, s+h2, inv_brho)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k3 = kick_vector_normalized(x3, px3, y3, py3, z3, pz3, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  x4 = x+h*k3[1]; px4 = px+h*k3[2]; y4 = y+h*k3[3]
  py4 = py+h*k3[4]; z4 = z+h*k3[5]; pz4 = pz+h*k3[6]
  Ex1, Ey1, Ez1, Bx1, By1, Bz1 = eval_mp_normalized(mp_field, x4, y4, z4, pz4, s+h)
  Ex2, Ey2, Ez2, Bx2, By2, Bz2 = eval_fm_normalized(fm_field, x4, y4, z4, pz4, s+h, inv_brho)
  Ex = Ex1+Ex2; Ey = Ey1+Ey2; Ez = Ez1+Ez2
  Bx = Bx1+Bx2; By = By1+By2; Bz = Bz1+Bz2
  k4 = kick_vector_normalized(x4, px4, y4, py4, z4, pz4, s+h, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  h6 = h / 6
  v[i, XI]  = vifelse(alive, x  + h6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]), v[i, XI])
  v[i, PXI] = vifelse(alive, px + h6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]), v[i, PXI])
  v[i, YI]  = vifelse(alive, y  + h6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]), v[i, YI])
  v[i, PYI] = vifelse(alive, py + h6*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]), v[i, PYI])
  v[i, ZI]  = vifelse(alive, z  + h6*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5]), v[i, ZI])
  v[i, PZI] = vifelse(alive, pz + h6*(k1[6] + 2*k2[6] + 2*k3[6] + k4[6]), v[i, PZI])
end

@makekernel function rk4_kernel_super_norm!(i, coords::Coords, beta_0, tilde_m,
                                p0c, mc2, s_span, ds_step, mp_field, fm_field, inv_brho)
  s_start = s_span[1]; s_end = s_span[2]; s = s_start
  v = coords.v
  total_distance = s_end - s_start
  n_steps = ceil(Int, total_distance / ds_step)
  for step in 1:n_steps
    remaining = s_end - s; h = min(ds_step, remaining)
    rel_p = 1 + v[i, PZI]; inv_rel_p = 1 / rel_p
    vt2 = (v[i, PXI] * inv_rel_p)^2 + (v[i, PYI] * inv_rel_p)^2
    alive = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse((vt2 >= 1) & alive, STATE_LOST_PZ, coords.state[i])
    rk4_step_super_norm!(coords, i, s, h, mp_field, fm_field, tilde_m, beta_0, p0c, mc2, inv_brho)
    s += h
  end
end
end # SuperpositionNormalized

# ============================================================
# Normalized-units kernels (proposed by reviewers)
# ============================================================
module NormalizedKernel
using BeamTracking: @makekernel, Coords
using BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using BeamTracking: C_LIGHT, vifelse, normalized_field, fieldmap_em_field, FieldMap, RectGrid3D
using BeamTracking.RungeKuttaTracking: MultipoleSource, FieldMapSource
using StaticArrays

# --- Normalized field evaluation ---

# Multipole: skip p_over_q_ref multiply (saves 3 multiplies per eval)
@inline function eval_em_field_normalized(f::MultipoleSource, x, y, z, pz, s)
  mm = f.mm; kn = f.kn; ks = f.ks
  return _multipole_normalized(x, y, z, s, mm, kn, ks)
end

@inline function _multipole_normalized(x, y, z, s, mm::SVector{0}, kn, ks)
  return (zero(x), zero(x), zero(x), zero(x), zero(x), zero(x))
end

@inline function _multipole_normalized(x, y, z, s, mm::SVector{N}, kn, ks) where N
  bx, by = normalized_field(mm, kn, ks, x, y, 0)
  is_solenoid = (mm[1] == 0)
  bz = vifelse(is_solenoid, kn[1], zero(x))
  return (zero(x), zero(x), zero(x), bx, by, bz)  # no * p_over_q_ref
end

# Field map: divide by p_over_q_ref after interpolation (adds 6 multiplies per eval)
@inline function eval_em_field_normalized(f::FieldMapSource, x, y, z, pz, s, inv_brho)
  Ex, Ey, Ez, Bx, By, Bz = fieldmap_em_field(x, y, z, pz, s, f.fieldmap, f.z_offset)
  return (Ex * inv_brho, Ey * inv_brho, Ez * inv_brho,
          Bx * inv_brho, By * inv_brho, Bz * inv_brho)
end

# For multipole, inv_brho is unused but we accept it for uniform dispatch
@inline function eval_em_field_normalized(f::MultipoleSource, x, y, z, pz, s, inv_brho)
  return eval_em_field_normalized(f, x, y, z, pz, s)
end

@inline get_g_bend(f::MultipoleSource) = f.g_bend
@inline get_g_bend(f::FieldMapSource) = f.g_bend

# --- kick_vector with simplified EOM using normalized fields ---
# With en = E/Bρ, bn = B/Bρ, the factors charge/p0 = 1/Bρ cancel:
#   dpx/ds = (en_x + vy*bn_z - vz*bn_y) * dt/ds  (no charge, no /p0)
#   (v×B)·v = 0, so F·v = qE·v only, and dp/ds simplifies similarly

@inline function kick_vector_normalized(x, px, y, py, z, pz, s, en_x, en_y, en_z, bn_x, bn_y, bn_z,
                tilde_m, beta_0, g_bend, p0c, mc2)
  rel_p = 1 + pz
  vt_x = px / rel_p
  vt_y = py / rel_p
  vt2 = vt_x^2 + vt_y^2

  vt2_1 = one(vt2)
  good_momenta = (vt2 < vt2_1)
  vt2_safe = vifelse(good_momenta, vt2, zero(vt2))

  rel_p2 = rel_p^2
  inv_gamma_v = sqrt(rel_p2 + tilde_m^2)
  beta = rel_p / inv_gamma_v
  inv_beta_c = 1 / (beta * C_LIGHT)

  rel_dir = 1
  vz_norm = sqrt(1 - vt2_safe) * rel_dir
  vx = beta * C_LIGHT * vt_x
  vy = beta * C_LIGHT * vt_y
  vz = beta * C_LIGHT * vz_norm

  # Normalized force: no charge, no /p0 — they cancel with Bρ
  # F_norm = en + v × bn  (force per Bρ, already divided out)
  Fn_x = en_x + (vy*bn_z - vz*bn_y)
  Fn_y = en_y + (vz*bn_x - vx*bn_z)

  dh_bend = x * g_bend
  abs_vz = abs(vz)
  abs_vz_safe = vifelse(good_momenta, abs_vz, one(abs_vz))
  dt_ds = rel_dir * (1 + dh_bend) / abs_vz_safe

  pz_p0 = rel_p * rel_dir * abs_vz * inv_beta_c

  # Energy: (v×B)·v = 0, so F·v = qE·v. Normalized: en·v / (beta*c) = dpz/ds
  en_dot_v = en_x*vx + en_y*vy + en_z*vz
  dpz_ds_val = en_dot_v * dt_ds * inv_beta_c

  # dbeta_ds needs physical dp/ds = dpz_ds * p0 = dpz_ds * p0c/c
  e_tot = p0c * rel_p / beta
  dbeta_ds = mc2^2 * dpz_ds_val * p0c / e_tot^3

  dx_ds = vx * dt_ds
  dy_ds = vy * dt_ds

  # dpx/ds = Fn_x * dt/ds  (charge/p0 already cancelled)
  dpx_ds = Fn_x * dt_ds + g_bend * pz_p0
  dpy_ds = Fn_y * dt_ds

  sqrt_1mvt2 = sqrt(1 - vt2_safe)
  dz_ds = rel_dir * (beta / beta_0 - 1) + rel_dir * (sqrt_1mvt2 - 1 - dh_bend) / sqrt_1mvt2 + dbeta_ds * z / beta

  zero_deriv = zero(dx_ds)
  return SVector(
    vifelse(good_momenta, dx_ds, zero_deriv),
    vifelse(good_momenta, dpx_ds, zero_deriv),
    vifelse(good_momenta, dy_ds, zero_deriv),
    vifelse(good_momenta, dpy_ds, zero_deriv),
    vifelse(good_momenta, dz_ds, zero_deriv),
    vifelse(good_momenta, dpz_ds_val, zero_deriv)
  )
end

# --- Normalized RK4 step ---

@inline function rk4_step_normalized!(coords, i, s, h, field, tilde_m, beta_0, p0c, mc2, inv_brho)
  g_bend = get_g_bend(field)
  alive = (coords.state[i] == STATE_ALIVE)
  v = coords.v
  x  = v[i, XI]; px = v[i, PXI]
  y  = v[i, YI]; py = v[i, PYI]
  z  = v[i, ZI]; pz = v[i, PZI]

  Ex, Ey, Ez, Bx, By, Bz = eval_em_field_normalized(field, x, y, z, pz, s, inv_brho)
  k1 = kick_vector_normalized(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  h2 = h / 2
  x2 = x+h2*k1[1]; px2 = px+h2*k1[2]; y2 = y+h2*k1[3]
  py2 = py+h2*k1[4]; z2 = z+h2*k1[5]; pz2 = pz+h2*k1[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field_normalized(field, x2, y2, z2, pz2, s+h2, inv_brho)
  k2 = kick_vector_normalized(x2, px2, y2, py2, z2, pz2, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  x3 = x+h2*k2[1]; px3 = px+h2*k2[2]; y3 = y+h2*k2[3]
  py3 = py+h2*k2[4]; z3 = z+h2*k2[5]; pz3 = pz+h2*k2[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field_normalized(field, x3, y3, z3, pz3, s+h2, inv_brho)
  k3 = kick_vector_normalized(x3, px3, y3, py3, z3, pz3, s+h2, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  x4 = x+h*k3[1]; px4 = px+h*k3[2]; y4 = y+h*k3[3]
  py4 = py+h*k3[4]; z4 = z+h*k3[5]; pz4 = pz+h*k3[6]
  Ex, Ey, Ez, Bx, By, Bz = eval_em_field_normalized(field, x4, y4, z4, pz4, s+h, inv_brho)
  k4 = kick_vector_normalized(x4, px4, y4, py4, z4, pz4, s+h, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, g_bend, p0c, mc2)

  h6 = h / 6
  v[i, XI]  = vifelse(alive, x  + h6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]), v[i, XI])
  v[i, PXI] = vifelse(alive, px + h6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]), v[i, PXI])
  v[i, YI]  = vifelse(alive, y  + h6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]), v[i, YI])
  v[i, PYI] = vifelse(alive, py + h6*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]), v[i, PYI])
  v[i, ZI]  = vifelse(alive, z  + h6*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5]), v[i, ZI])
  v[i, PZI] = vifelse(alive, pz + h6*(k1[6] + 2*k2[6] + 2*k3[6] + k4[6]), v[i, PZI])
end

# --- Normalized RK4 kernel ---

@makekernel function rk4_kernel_normalized!(i, coords::Coords, beta_0, tilde_m,
                                p0c, mc2, s_span, ds_step, field, inv_brho)
  s_start = s_span[1]; s_end = s_span[2]; s = s_start
  v = coords.v
  total_distance = s_end - s_start
  n_steps = ceil(Int, total_distance / ds_step)
  for step in 1:n_steps
    remaining = s_end - s; h = min(ds_step, remaining)
    rel_p = 1 + v[i, PZI]; inv_rel_p = 1 / rel_p
    vt2 = (v[i, PXI] * inv_rel_p)^2 + (v[i, PYI] * inv_rel_p)^2
    alive = (coords.state[i] == STATE_ALIVE)
    coords.state[i] = vifelse((vt2 >= 1) & alive, STATE_LOST_PZ, coords.state[i])
    rk4_step_normalized!(coords, i, s, h, field, tilde_m, beta_0, p0c, mc2, inv_brho)
    s += h
  end
end

end # NormalizedKernel

# ============================================================
# Setup helpers
# ============================================================

function setup_particle_params(pc=1e9)
  species = Species("electron")
  mc2 = massof(species)
  p_over_q_ref = pc_to_R(species, pc)
  beta_gamma_0 = R_to_beta_gamma(species, p_over_q_ref)
  tilde_m = 1 / beta_gamma_0
  beta_0 = beta_gamma_0 / sqrt(1 + beta_gamma_0^2)
  charge = chargeof(species)
  p0c = R_to_pc(species, p_over_q_ref)
  return species, p_over_q_ref, beta_0, tilde_m, charge, p0c, mc2
end

function make_bunch(n_particles, species, p_over_q_ref)
  Bunch(randn(n_particles, 6) * 0.001, p_over_q_ref=p_over_q_ref, species=species)
end

function make_multipole_source(p_over_q_ref)
  mm = SVector(1, 2)        # dipole + quadrupole
  kn = SVector(0.0, 0.5)    # normalized strengths
  ks = SVector(0.0, 0.0)
  g_bend = 0.0
  MultipoleSource(mm, kn, ks, p_over_q_ref, g_bend)
end

function make_fieldmap_source(p_over_q_ref)
  # 3D rectangular grid: quadrupole-like B field in physical units (Tesla)
  nx, ny, nz = 20, 20, 50
  dx, dy, dz = 0.001, 0.001, 0.02  # 2cm x 2cm x 1m
  grid = RectGrid3D(
    (-nx/2*dx, -ny/2*dy, 0.0),
    (dx, dy, dz),
    (0, 0, 0),
    (nx, ny, nz)
  )
  k_quad = 0.5 * p_over_q_ref  # normalized strength 0.5/m -> physical Tesla/m
  B = zeros(3, nx, ny, nz)
  for iz in 1:nz, iy in 1:ny, ix in 1:nx
    x = grid.gridOriginOffset[1] + (ix-1)*dx
    y = grid.gridOriginOffset[2] + (iy-1)*dy
    B[1, ix, iy, iz] = k_quad * y   # Bx = k*y
    B[2, ix, iy, iz] = k_quad * x   # By = k*x
  end
  fm = FieldMap(grid; B=B)
  FieldMapSource(fm, 0.0, 0.0)
end

function make_chain(kernel, params)
  kc = KernelChain(Val{1}())
  kc = BeamTracking.push(kc, KernelCall(kernel, params))
  return kc
end

function reset_coords!(coords, v_init, state_init)
  coords.v .= v_init
  coords.state .= state_init
end

# ============================================================
# Benchmark runner
# ============================================================

function report(label, b)
  t_min = minimum(b).time / 1e3
  t_med = median(b).time / 1e3
  t_std = std(b.times) / 1e3
  n = length(b.times)
  println("    $label: min=$(round(t_min, digits=2)) med=$(round(t_med, digits=2)) ±$(round(t_std, digits=2)) μs (n=$n)")
  return median(b).time
end

function run_benchmark_cpu(label, coords, kc, v_init, state_init; use_simd=true)
  reset_coords!(coords, v_init, state_init)
  launch!(coords, kc; use_KA=false, use_explicit_SIMD=use_simd)

  b = @benchmark begin
    launch!($coords, $kc; use_KA=false, use_explicit_SIMD=$use_simd)
  end setup=(reset_coords!($coords, $v_init, $state_init))

  return report(label, b)
end

function run_benchmark_gpu(label, coords_gpu, kc, v_init_gpu, state_init_gpu)
  coords_gpu.v .= v_init_gpu
  coords_gpu.state .= state_init_gpu
  launch!(coords_gpu, kc; use_KA=true, use_explicit_SIMD=false)

  b = @benchmark begin
    $coords_gpu.v .= $v_init_gpu
    $coords_gpu.state .= $state_init_gpu
    launch!($coords_gpu, $kc; use_KA=true, use_explicit_SIMD=false)
  end

  return report(label, b)
end

function print_diff(t_phys, t_norm)
  pct = (t_phys - t_norm) / t_phys * 100
  println("    -> Δ = $(round(pct, digits=2))% (positive = physical slower)\n")
end

# ============================================================
# Main
# ============================================================

function main()
  species, p_over_q_ref, beta_0, tilde_m, charge, p0c, mc2 = setup_particle_params()

  mp_source = make_multipole_source(p_over_q_ref)
  fm_source = make_fieldmap_source(p_over_q_ref)

  s_span = (0.0, 1.0)
  ds_step = 0.01

  # Physical (current): uses RungeKuttaTracking.rk4_kernel! with field source structs
  params_mp_phys = (beta_0, tilde_m, charge, p0c, mc2, s_span, ds_step, mp_source)
  params_fm_phys = (beta_0, tilde_m, charge, p0c, mc2, s_span, ds_step, fm_source)

  # Normalized (proposed): simplified EOM, no charge/p0 in force terms
  inv_brho = 1 / p_over_q_ref
  params_mp_norm = (beta_0, tilde_m, p0c, mc2, s_span, ds_step, mp_source, inv_brho)
  params_fm_norm = (beta_0, tilde_m, p0c, mc2, s_span, ds_step, fm_source, inv_brho)

  kc_mp_phys = make_chain(RungeKuttaTracking.rk4_kernel!, params_mp_phys)
  kc_mp_norm = make_chain(NormalizedKernel.rk4_kernel_normalized!, params_mp_norm)
  kc_fm_phys = make_chain(RungeKuttaTracking.rk4_kernel!, params_fm_phys)
  kc_fm_norm = make_chain(NormalizedKernel.rk4_kernel_normalized!, params_fm_norm)

  for n_particles in [1_000, 10_000, 100_000]
    bunch = make_bunch(n_particles, species, p_over_q_ref)
    coords = bunch.coords
    v_init = copy(coords.v)
    state_init = copy(coords.state)

    println("="^70)
    println("N = $n_particles particles | 100 RK4 steps | dipole+quad")
    println("="^70)

    # --- Correctness check ---
    println("\n  Correctness check:")
    for (label, kc_phys, kc_norm) in [
        ("Multipole", kc_mp_phys, kc_mp_norm),
        ("Field map", kc_fm_phys, kc_fm_norm)]
      reset_coords!(coords, v_init, state_init)
      launch!(coords, kc_phys; use_KA=false, use_explicit_SIMD=false)
      v_phys = copy(coords.v)

      reset_coords!(coords, v_init, state_init)
      launch!(coords, kc_norm; use_KA=false, use_explicit_SIMD=false)
      v_norm = copy(coords.v)

      max_diff = maximum(abs.(v_phys .- v_norm))
      rel_diff = maximum(abs.(v_phys .- v_norm) ./ (abs.(v_phys) .+ 1e-30))
      println("    $label: max abs diff = $max_diff, max rel diff = $rel_diff")
      reset_coords!(coords, v_init, state_init)
    end

    # --- CPU ---
    println("\n  CPU — Analytical multipole:")
    t1 = run_benchmark_cpu("Physical  (current) ", coords, kc_mp_phys, v_init, state_init)
    t2 = run_benchmark_cpu("Normalized (proposed)", coords, kc_mp_norm, v_init, state_init)
    print_diff(t1, t2)

    println("  CPU — Field map (3D rect grid, no SIMD):")
    t1 = run_benchmark_cpu("Physical  (current) ", coords, kc_fm_phys, v_init, state_init; use_simd=false)
    t2 = run_benchmark_cpu("Normalized (proposed)", coords, kc_fm_norm, v_init, state_init; use_simd=false)
    print_diff(t1, t2)

    # --- GPU ---
    if _has_cuda
      v_gpu = CuArray(v_init)
      state_gpu = CuArray(state_init)
      coords_gpu = Coords(state_gpu, v_gpu, nothing)
      v_init_gpu = CuArray(v_init)
      state_init_gpu = CuArray(state_init)

      println("  GPU (CUDA) — Analytical multipole:")
      t1 = run_benchmark_gpu("Physical  (current) ", coords_gpu, kc_mp_phys, v_init_gpu, state_init_gpu)
      t2 = run_benchmark_gpu("Normalized (proposed)", coords_gpu, kc_mp_norm, v_init_gpu, state_init_gpu)
      print_diff(t1, t2)

      # Field map on GPU: adapt the FieldMap, then rebuild FieldMapSource
      fm_gpu = Adapt.adapt(CuArray, fm_source.fieldmap)
      fm_source_gpu = FieldMapSource(fm_gpu, fm_source.z_offset, fm_source.g_bend)
      params_fm_phys_gpu = (beta_0, tilde_m, charge, p0c, mc2, s_span, ds_step, fm_source_gpu)
      params_fm_norm_gpu = (beta_0, tilde_m, p0c, mc2, s_span, ds_step, fm_source_gpu, inv_brho)
      kc_fm_phys_gpu = make_chain(RungeKuttaTracking.rk4_kernel!, params_fm_phys_gpu)
      kc_fm_norm_gpu = make_chain(NormalizedKernel.rk4_kernel_normalized!, params_fm_norm_gpu)

      println("  GPU (CUDA) — Field map (3D rect grid):")
      t1 = run_benchmark_gpu("Physical  (current) ", coords_gpu, kc_fm_phys_gpu, v_init_gpu, state_init_gpu)
      t2 = run_benchmark_gpu("Normalized (proposed)", coords_gpu, kc_fm_norm_gpu, v_init_gpu, state_init_gpu)
      print_diff(t1, t2)
    else
      println("  GPU: CUDA not available, skipping\n")
    end

    # --- Superposition: multipole + field map ---
    params_super_phys = (beta_0, tilde_m, charge, p0c, mc2, s_span, ds_step, mp_source, fm_source)
    params_super_norm = (beta_0, tilde_m, p0c, mc2, s_span, ds_step, mp_source, fm_source, inv_brho)
    kc_super_phys = make_chain(SuperpositionPhysical.rk4_kernel_super!, params_super_phys)
    kc_super_norm = make_chain(SuperpositionNormalized.rk4_kernel_super_norm!, params_super_norm)

    # Correctness
    reset_coords!(coords, v_init, state_init)
    launch!(coords, kc_super_phys; use_KA=false, use_explicit_SIMD=false)
    v_phys = copy(coords.v)
    reset_coords!(coords, v_init, state_init)
    launch!(coords, kc_super_norm; use_KA=false, use_explicit_SIMD=false)
    v_norm = copy(coords.v)
    max_diff = maximum(abs.(v_phys .- v_norm))
    rel_diff = maximum(abs.(v_phys .- v_norm) ./ (abs.(v_phys) .+ 1e-30))
    println("  Superposition correctness: max abs = $max_diff, max rel = $rel_diff")
    reset_coords!(coords, v_init, state_init)

    println("\n  CPU — Superposition (multipole + field map, no SIMD):")
    t1 = run_benchmark_cpu("Physical  (current) ", coords, kc_super_phys, v_init, state_init; use_simd=false)
    t2 = run_benchmark_cpu("Normalized (proposed)", coords, kc_super_norm, v_init, state_init; use_simd=false)
    print_diff(t1, t2)

    if _has_cuda
      fm_gpu_s = Adapt.adapt(CuArray, fm_source.fieldmap)
      fm_source_gpu_s = FieldMapSource(fm_gpu_s, fm_source.z_offset, fm_source.g_bend)
      params_super_phys_gpu = (beta_0, tilde_m, charge, p0c, mc2, s_span, ds_step, mp_source, fm_source_gpu_s)
      params_super_norm_gpu = (beta_0, tilde_m, p0c, mc2, s_span, ds_step, mp_source, fm_source_gpu_s, inv_brho)
      kc_super_phys_gpu = make_chain(SuperpositionPhysical.rk4_kernel_super!, params_super_phys_gpu)
      kc_super_norm_gpu = make_chain(SuperpositionNormalized.rk4_kernel_super_norm!, params_super_norm_gpu)

      println("  GPU (CUDA) — Superposition (multipole + field map):")
      t1 = run_benchmark_gpu("Physical  (current) ", coords_gpu, kc_super_phys_gpu, v_init_gpu, state_init_gpu)
      t2 = run_benchmark_gpu("Normalized (proposed)", coords_gpu, kc_super_norm_gpu, v_init_gpu, state_init_gpu)
      print_diff(t1, t2)
    end
    println()
  end
end

main()
