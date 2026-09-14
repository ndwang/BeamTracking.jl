"""
  RungeKuttaTracking

Module implementing particle tracking through electromagnetic field sources
using a fourth-order Runge-Kutta method.
"""
module RungeKuttaTracking
using ..BeamTracking, ..StaticArrays
using ..BeamTracking: @makekernel, Coords
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, STATE_ALIVE, STATE_LOST_PZ
using ..BeamTracking: C_LIGHT, EMField, vifelse

"""
  _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
        tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

Calculate the derivative vector du/ds for relativistic particle tracking.
Returns an SVector{6} containing [dx/ds, dpx/ds, dy/ds, dpy/ds, dz/ds, dpz/ds].

Uses branchless operations for GPU/SIMD compatibility. For unphysical momenta,
returns zero derivatives (caller should mark particle as lost).

# Arguments
- `x, px, y, py, z, pz`: State vector components
- `s`: Arc length position
- `Ex, Ey, Ez`: Physical electric field (V/m)
- `Bx, By, Bz`: Physical magnetic field (T)
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gx`, `gy`: Horizontal and vertical reference curvature components
- `electric_scale`, `magnetic_scale`: q/(p₀c) and q/p₀ conversion factors
"""
@inline function _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
  # Relative momentum
  rel_p = 1 + pz

  # Transverse velocity components (normalized)
  inv_rel_p = inv(rel_p)
  vt_x = px * inv_rel_p
  vt_y = py * inv_rel_p
  vt2 = vt_x^2 + vt_y^2

  # Check for unphysical momenta (branchless)
  vt2_1 = one(vt2)
  good_momenta = (vt2 < vt2_1)
  vt2_safe = vifelse(good_momenta, vt2, zero(vt2))

  # Particle beta and velocity
  rel_p2 = rel_p^2
  inv_gamma_v = sqrt(rel_p2 + tilde_m^2)
  beta = rel_p / inv_gamma_v
  
  # Unit direction and path-length factor; field conversion factors are
  # computed once per particle, outside the RK stages.
  uz = sqrt(1 - vt2_safe)
  dh_bend = x * gx + y * gy
  inv_uz = inv(uz)
  path_factor = (1 + dh_bend) * inv_uz
  inv_beta = inv_gamma_v * inv_rel_p

  dx_ds = vt_x * path_factor
  dy_ds = vt_y * path_factor
  dpx_ds = (Ex * electric_scale * inv_beta + (vt_y * Bz - uz * By) * magnetic_scale) * path_factor + gx * rel_p * uz
  dpy_ds = (Ey * electric_scale * inv_beta + (uz * Bx - vt_x * Bz) * magnetic_scale) * path_factor + gy * rel_p * uz

  # Only E changes |p|. dβ/ds = (m c/p₀)² / (E/(p₀ c))³ * d(|p|/p₀)/ds.
  dpz_ds = (Ex * vt_x + Ey * vt_y + Ez * uz) * electric_scale * path_factor * inv_beta
  # dβ/ds * z/β = (m̃² / (rel_p² + m̃²)) * dpz/ds * z/rel_p.
  dz_accel = (tilde_m^2 / (rel_p2 + tilde_m^2)) * dpz_ds * z * inv_rel_p
  dz_ds = beta / beta_0 - path_factor + dz_accel

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

@inline function _kick_vector(x, px, y, py, z, pz, s, field::EMField,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
  Ex, Ey, Ez = field.E
  Bx, By, Bz = field.B
  return _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                     tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
end

"""
    kick_vector(x, px, y, py, z, pz, s, field::EMField,
                charge, tilde_m, beta_0, gx, gy, p0c, mc2)

Evaluate the simplified RK equations with physical electric (V/m) and magnetic
(T) fields. `charge` is in units of e and `p0c` is in eV. The original argument
list is preserved; `mc2` is redundant with `tilde_m` and `p0c` and is unused.
The component overload accepts `Ex, Ey, Ez, Bx, By, Bz` in place of `field`.
"""
@inline function kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                             charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  electric_scale = charge / p0c
  magnetic_scale = electric_scale * C_LIGHT
  return _kick_vector(x, px, y, py, z, pz, s, Ex, Ey, Ez, Bx, By, Bz,
                      tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
end

@inline function kick_vector(x, px, y, py, z, pz, s, field::EMField,
                             charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  return kick_vector(x, px, y, py, z, pz, s, field.E..., field.B...,
                     charge, tilde_m, beta_0, gx, gy, p0c, mc2)
end

"""
  _rk4_step!(coords, i, s, h, source, tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

Perform a single RK4 step for particle i, updating coordinates in-place.
Only updates state if particle is alive.

# Arguments
- `coords`: Coordinates structure
- `i`: Particle index
- `s`: Current arc length
- `h`: Step size
- `source`: Concrete callable field source
- `tilde_m`: Normalized mass mc²/(p₀c)
- `beta_0`: Reference velocity β₀ = v₀/c
- `gx`, `gy`: Horizontal and vertical reference curvature components
- `electric_scale`, `magnetic_scale`: q/(p₀c) and q/p₀ conversion factors
"""
@inline function _rk4_step!(coords, i, s, h, source, tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
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
  field = source(x, y, z, s)
  k1 = _kick_vector(x, px, y, py, z, pz, s, field,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

  # k2 = f(u + h/2 * k1, s + h/2)
  h2 = h / 2
  x2 = x + h2 * k1[1]
  px2 = px + h2 * k1[2]
  y2 = y + h2 * k1[3]
  py2 = py + h2 * k1[4]
  z2 = z + h2 * k1[5]
  pz2 = pz + h2 * k1[6]
  field = source(x2, y2, z2, s + h2)
  k2 = _kick_vector(x2, px2, y2, py2, z2, pz2, s + h2, field,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

  # k3 = f(u + h/2 * k2, s + h/2)
  x3 = x + h2 * k2[1]
  px3 = px + h2 * k2[2]
  y3 = y + h2 * k2[3]
  py3 = py + h2 * k2[4]
  z3 = z + h2 * k2[5]
  pz3 = pz + h2 * k2[6]
  field = source(x3, y3, z3, s + h2)
  k3 = _kick_vector(x3, px3, y3, py3, z3, pz3, s + h2, field,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

  # k4 = f(u + h * k3, s + h)
  x4 = x + h * k3[1]
  px4 = px + h * k3[2]
  y4 = y + h * k3[3]
  py4 = py + h * k3[4]
  z4 = z + h * k3[5]
  pz4 = pz + h * k3[6]
  field = source(x4, y4, z4, s + h)
  k4 = _kick_vector(x4, px4, y4, py4, z4, pz4, s + h, field,
                tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)

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
    rk4_step!(coords, i, s, h, source, charge, tilde_m, beta_0, gx, gy, p0c, mc2)

Advance one RK4 step using a physical-unit field source. The full tracking
kernel reuses the field conversion factors across all steps for each particle.
"""
@inline function rk4_step!(coords, i, s, h, source, charge, tilde_m, beta_0, gx, gy, p0c, mc2)
  electric_scale = charge / p0c
  return _rk4_step!(coords, i, s, h, source, tilde_m, beta_0, gx, gy,
                   electric_scale, electric_scale * C_LIGHT)
end

"""
  rk4_kernel!(i, coords, beta_0, tilde_m, charge, p0c, mc2,
              L, ds_step, n_steps, gx, gy, source)

Kernelized RK4 tracking through a concrete electromagnetic field source.
Compatible with @makekernel and the package's kernel architecture.
"""
@makekernel function rk4_kernel!(i, coords::Coords, beta_0, tilde_m, charge, p0c, mc2,
                                L, ds_step, n_steps,
                                gx, gy, source)
  s = zero(L)
  electric_scale = charge / p0c
  magnetic_scale = electric_scale * C_LIGHT

  v = coords.v
  
  for step in 1:n_steps
    # Check if particle is lost
    rel_p = 1 + v[i, PZI]
    inv_rel_p = 1 / rel_p
    vt2 = (v[i, PXI] * inv_rel_p)^2 + (v[i, PYI] * inv_rel_p)^2
    alive = (coords.state[i] == STATE_ALIVE)
    # Mark particle as lost
    coords.state[i] = vifelse((vt2 >= 1) & alive, STATE_LOST_PZ, coords.state[i])

    # Perform RK4 step (check for alive status is now inside rk4_step!)
    _rk4_step!(coords, i, s, ds_step, source, tilde_m, beta_0, gx, gy, electric_scale, magnetic_scale)
    s += ds_step

    # The common path performs the final callback after exit processing.
    if step != n_steps
      BeamTracking.execute_callbacks(i, coords, s, s / (beta_0 * C_LIGHT))
    end
  end
end

end
