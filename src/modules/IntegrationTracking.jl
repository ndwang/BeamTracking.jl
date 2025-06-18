#=

Tracking using symplectic integration.

=#

@Base.kwdef struct Integration
  order::Int32 = 2
  n_steps::Int32 = 1
end

module IntegrationTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays, ..KernelAbstractions
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, @makekernel, BunchView
const TRACKING_METHOD = Integration

@makekernel function order_two_integrator!(i, b::BunchView, ker, params, ds_step, n_steps, L)
  for _ in 1:n_steps
    ker(i, b, params..., ds_step)
  end
end


@makekernel function order_four_integrator!(i, b::BunchView, ker, params, ds_step, n_steps, L)
  w0 = 1.3512071919596577718181151794851757586002349853515625*ds_step
  w1 = -1.7024143839193153215916254339390434324741363525390625*ds_step
  for _ in 1:n_steps
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
  end
end


@makekernel function order_six_integrator!(i, b, ker, params, ds_step, n_steps, L)
  w0 = 1.3151863206857402*ds_step
  w1 = -1.17767998417887*ds_step
  w2 = 0.235573213359*ds_step
  w3 = 0.784513610477*ds_step
  for _ in 1:n_steps
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
  end
end


@makekernel function order_eight_integrator!(i, b, ker, params, ds_step, n_steps, L)
  w0 = -1.7808286265894515*ds_step
  w1 = -1.61582374150097*ds_step
  w2 = -2.44699182370524*ds_step
  w3 = -0.716989419708120e-2*ds_step
  w4 = 2.44002732616735*ds_step
  w5 = 0.157739928123617*ds_step
  w6 = 1.82020630970714*ds_step
  w7 = 1.04242620869991*ds_step
  for _ in 1:n_steps
    ker(i, b, params..., w7)
    ker(i, b, params..., w6)
    ker(i, b, params..., w5)
    ker(i, b, params..., w4)
    ker(i, b, params..., w3)
    ker(i, b, params..., w2)
    ker(i, b, params..., w1)
    ker(i, b, params..., w0)
    ker(i, b, params..., w1)
    ker(i, b, params..., w2)
    ker(i, b, params..., w3)
    ker(i, b, params..., w4)
    ker(i, b, params..., w5)
    ker(i, b, params..., w6)
    ker(i, b, params..., w7)
  end
end

end