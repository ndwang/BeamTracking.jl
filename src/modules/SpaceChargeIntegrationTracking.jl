struct SpaceChargeIntegration
  order::Int
  num_steps::Int
  num_sc_steps::Int
  ds_step::Float64

  function SpaceChargeIntegration(; order::Int=2, num_steps::Int=-1, num_sc_steps::Int=-1, ds_step::Float64=-1.0)
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
      return new(_order, _num_steps, num_sc_steps, _ds_step)
  end
end

module SpaceChargeIntegrationTracking
using ..GTPSA, ..BeamTracking, ..KernelAbstractions, ..SpaceCharge
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, Q0, QX, QY, QZ, @makekernel, Coords
const TRACKING_METHOD = SpaceChargeIntegration


function sc_calc(scp, bunch)
  # Extract only alive particles
  alive_idx = findall(==(State.Alive), bunch.coords.state)
  particles_x = bunch.coords.v[alive_idx, XI]
  particles_y = bunch.coords.v[alive_idx, YI]
  particles_z = bunch.coords.v[alive_idx, ZI]
  particles_q = scp.total_charge * Float64.(bunch.coords.state .== State.Alive) / length(bunch.coords.state)
  scp.mesh.gamma  = BeamTracking.R_to_gamma(bunch.species, bunch.R_ref * (1 + sum(bunch.coords.v[alive_idx, PZI])/length(alive_idx)))

  deposit!(scp.mesh, particles_x, particles_y, particles_z, particles_q)
  solve!(scp.mesh)
end

@makekernel fastgtpsa=true function interpolate_field(i, coords::Coords, scp)
  v = coords.v
  efield = scp.mesh.efield
  min_bounds = scp.mesh.min_bounds
  delta = scp.mesh.delta

  norm_x = (v[i, XI] - min_bounds[1]) / delta[1]
  norm_y = (v[i, YI] - min_bounds[2]) / delta[2]
  norm_z = (v[i, ZI] - min_bounds[3]) / delta[3]

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

  scp.efield_scratch[i, 1] = (
    efield[ix + 1, iy + 1, iz + 1, 1] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 1] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 1] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 1] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 1] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 1] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 1] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 1] * w_111
  )
  scp.efield_scratch[i, 2] = (
    efield[ix + 1, iy + 1, iz + 1, 2] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 2] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 2] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 2] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 2] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 2] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 2] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 2] * w_111
  )
  scp.efield_scratch[i, 3] = (
    efield[ix + 1, iy + 1, iz + 1, 3] * w_000 +
    efield[ix + 2, iy + 1, iz + 1, 3] * w_100 +
    efield[ix + 1, iy + 2, iz + 1, 3] * w_010 +
    efield[ix + 2, iy + 2, iz + 1, 3] * w_110 +
    efield[ix + 1, iy + 1, iz + 2, 3] * w_001 +
    efield[ix + 2, iy + 1, iz + 2, 3] * w_101 +
    efield[ix + 1, iy + 2, iz + 2, 3] * w_011 +
    efield[ix + 2, iy + 2, iz + 2, 3] * w_111
  )
end

@makekernel fastgtpsa=true function sc_kick!(i, coords::Coords, β_0, R, scp, L)
  v = coords.v
  Efield = scp.efield_scratch
  gamma2 = scp.mesh.gamma^2
  coeff = L / (β_0 * R * gamma2 * BeamTracking.C_LIGHT)

  Ps2 = (1 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2)
  coords.state[i] = ifelse(Ps2 <= 0 && coords.state[i] == State.Alive, State.Lost, coords.state[i])
  alive = ifelse(coords.state[i] == State.Alive, 1, 0) 

  Ps2_safe = max(0, Ps2)
  Ps = alive * (sqrt(Ps2_safe) + Efield[i,3] * coeff) + (1 - alive) * 0
  
  v[i, PXI] = alive * (v[i, PXI] + Efield[i,1] * coeff) + (1 - alive) * v[i, PXI]
  v[i, PYI] = alive * (v[i, PYI] + Efield[i,2] * coeff) + (1 - alive) * v[i, PYI]
  v[i, PZI] = alive * (sqrt(v[i, PXI]^2 + v[i, PYI]^2 + Ps^2) - 1) + (1 - alive) * v[i, PZI]
end

macro sc_wrap(kernel, params, beta_0, R, scp, L)
  quote
    function (i, coords, args...)
      sc_kick!(i, coords, $(esc(beta_0)), $(esc(R)), $(esc(scp)), $(esc(L))/2)
      $(esc(kernel))(i, coords, $(esc(params))..., args...)
      sc_kick!(i, coords, $(esc(beta_0)), $(esc(R)), $(esc(scp)), $(esc(L))/2)
    end
  end
end


end