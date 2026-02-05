module BeamTrackingCUDAExt
using CUDA: CuDeviceArray
import BeamTracking: gaussian_random

# CUDA.jl has horrifying error message with randn, but not rand
# Box-Muller transform let's us get around that
function gaussian_random(::CuDeviceArray, sigma1, sigma2)
  s, c = sincospi(2 * rand())
  t = sqrt(-2 * log(rand()))
  z0 = c*t*sigma1
  z1 = s*t*sigma2
  return z0, z1
end

end