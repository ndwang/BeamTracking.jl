using ..BeamTracking.IntegrationTracking.def_integrator_struct

struct SpaceCharge
  order::Int
  num_steps::Int
  num_sc_steps::Int
  ds_step::Float64

  function SpaceCharge(; order::Int=2, num_steps::Int=-1, num_sc_steps::Int=-1, ds_step::Float64=-1.0)
      _order = order
      _num_steps = num_steps
      _ds_step = ds_step
      if _order ∉ (2, 4, 6, 8)
          error("Symplectic integration only supports orders 2, 4, 6, and 8")
      elseif _num_steps == -1 && _ds_step == -1.0
          _num_steps = 1
      elseif _num_steps > 0 && _ds_step > 0
          error("Only one of num_steps or ds_step should be specified")
      elseif _num_steps < 1 && _ds_step <= 0
          error("Invalid step size")
      elseif _num_steps > 0
          _ds_step = -1.0
      elseif _ds_step > 0
          _num_steps = -1
      end
      return new(_order, _num_steps, _num_sc_steps, _ds_step)
  end
end

module SpaceChargeIntegration
using ..BeamTracking, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, @makekernel, Coords
const TRACKING_METHOD = SpaceCharge

@makekernel fastgtpsa=true function SC_kick!(i, coords::Coords, β_0, gamsqr_0, R, scp, L)
  efield = scp.efield

  norm_x = (coords.v[i, XI] - scp.min_bounds[1]) / scp.delta[1]
  norm_y = (coords.v[i, YI] - scp.min_bounds[2]) / scp.delta[2]
  norm_z = (coords.v[i, ZI] - scp.min_bounds[3]) / scp.delta[3]

  ix = floor(Int, norm_x)
  iy = floor(Int, norm_y)
  iz = floor(Int, norm_z)

  dx = norm_x - ix
  dy = norm_y - iy
  dz = norm_z - iz

  w_000 = (1 - dx) * (1 - dy) * (1 - dz)
  w_100 = dx * (1 - dy) * (1 - dz)
  w_010 = (1 - dx) * dy * (1 - dz)
  w_110 = dx * dy * (1 - dz)
  w_001 = (1 - dx) * (1 - dy) * dz
  w_101 = dx * (1 - dy) * dz
  w_011 = (1 - dx) * dy * dz
  w_111 = dx * dy * dz

  Ex = (
    efield[ix + 1, iy + 1, iz + 1, 1] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 1] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 1] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 1] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 1] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 1] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 1] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 1] * w_111
  )
  Ey = (
    efield[ix + 1, iy + 1, iz + 1, 2] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 2] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 2] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 2] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 2] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 2] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 2] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 2] * w_111
  )
  Ez = (
    efield[ix + 1, iy + 1, iz + 1, 3] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 3] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 3] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 3] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 3] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 3] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 3] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 3] * w_111
  )

  coeff = L * β_0 / (R * gamsqr_0)

  
  Ps2 = (1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2)
  coords.state[i] = ifelse(Ps2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i])
  alive = ifelse(coords.state[i]==State.Alive, 1, 0) 
  Ps = sqrt(Ps2 + (alive-1)*(Ps2-1))

  Ps              += Ez * coeff
  coords.v[i, PX] += Ex * coeff
  coords.v[i, PY] += Ey * coeff
  coords.v[i, PZ] = sqrt(coords.v[i, PX]^2 + coords.v[i, PY]^2 + Ps^2) - 1
end

@makekernel fastgtpsa=true function curved_drift_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, e1, e2, theta, g, w, w_inv, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  ExactTracking.exact_curved_drift!(i, coords, e1, e2, theta, g, w, w_inv, BeamTracking.anom(bunch.species), tilde_m, beta_0, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end

@makekernel fastgtpsa=true function drift_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  ExactTracking.exact_drift!(i, coords, beta_0, gamsqr_0, tilde_m, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end

@makekernel fastgtpsa=true function solenoid_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, Ksol, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  ExactTracking.exact_solenoid!(i, coords, Ksol, beta_0, gamsqr_0, tilde_m, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end

@makekernel fastgtpsa=true function bend_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, e1, e2, theta, g, Kn0, w, w_inv, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  ExactTracking.exact_bend!(i, coords, e1, e2, theta, g, Kn0, w, w_inv, tilde_m, beta_0, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end

@makekernel fastgtpsa=true function dkd_multipole_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, a, mm, kn, ks, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  IntegrationTracking.dkd_multipole!(i, coords, beta_0, gamsqr_0, tilde_m, a, mm, kn, ks, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end


@makekernel fastgtpsa=true function sks_multipole_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, a, Ksol, mm, kn, ks, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  IntegrationTracking.sks_multipole!(i, coords, beta_0, gamsqr_0, tilde_m, a, Ksol, mm, kn, ks, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end


@makekernel fastgtpsa=true function bkb_multipole_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, a, e1, e2, g, w, w_inv, k0, mm, kn, ks, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  IntegrationTracking.bkb_multipole!(i, coords, tilde_m, beta_0, a, e1, e2, g, w, w_inv, k0, mm, kn, ks, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end


@makekernel fastgtpsa=true function mkm_quadrupole_sc!(i, coords::Coords, tilde_m, gamsqr_0, beta_0, a, w, w_inv, k1, mm, kn, ks, R, scp, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
  IntegrationTracking.mkm_quadrupole!(i, coords, beta_0, gamsqr_0, tilde_m, a, w, w_inv, k1, mm, kn, ks, L)
  SC_kick!(i, coords::Coords, beta_0, gamsqr_0, R, scp, L/2)
end

end