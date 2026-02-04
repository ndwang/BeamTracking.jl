module BeamTrackingCUDAExt
using CUDA: CuArray
import BeamTracking: gaussian_random

function gaussian_random(::CuArray, sigma1, sigma2)
  s, c = sincospi(2 * rand())
  t = sqrt(-2 * log(rand()))
  z0 = c*t*sigma1
  z1 = s*t*sigma2
  return z0, z1
end

end