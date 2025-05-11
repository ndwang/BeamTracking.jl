
const REGISTER_SIZE = VectorizationBase.register_size()
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6


# Generic function to launch a kernel on the bunch coordinates matrix
# Matrix v should ALWAYS be in SoA whether for real or as a view via tranpose(v)

"""
    launch!(f!::F, v, v0, work, args...; simd_lane_width, multithread_threshold)

General purpose function to launch a kernel `f!`. The syntax for a kernel `f!` must 
ALWAYS be the following:

## Arguments
- `i`       -- Particle index
- `v`       -- Input/output matrix as an SoA or SoA view ALWAYS! (use transpose if AoS)
- `work`    -- A Vector of temporary vectors (columns of v) to run the kernel `f!`
- `args...` -- Any further arguments to run the kernel

## Keyword Arguments
- `simd_lane_width`       -- The number of SIMD lanes to use. Default is `REGISTER_SIZE/sizeof(eltype(A))`
- `multithread_threshold` -- Number of particles at which multithreading is used. Default is `1e6``
"""
@inline @generated function launch!(f!::F, v::V, args...) where {F<:Function,V}
  #=
  N_particle = size(v, 1)
  @simd for i in 1:N_particle
    @assert last(i) <= N_particle "Out of bounds!"
    f!(i, v, args...)
  end
=#
  
  backend = get_backend(similar(V,0,0))
  if backend isa CPU
    groupsize = floor(Int,REGISTER_SIZE/sizeof(eltype(V))) # Explicit SIMD
  else
    # Then this is CUDA block size
    groupsize = 256 
  end
  kernel! = tgeneric_kernel!(backend, groupsize)
  return quote
    N_particle = size(v, 1)
    $(kernel!)(f!, v, args; ndrange=N_particle)
    KernelAbstractions.synchronize($backend)
    return v
  end
  #=
  return quote

  end
  N_particle = size(v, 1)
  backend = get_backend(v)
  kernel! = tgeneric_kernel!(backend, groupsize)
  kernel!(f!, v, args; ndrange = N_particle)
  
  return v
  =#
  #=
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    else
      for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    end
    # Do the remainder
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, work, args...)
    end
  else
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    else
      @simd for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    end
  end
  return v
  =#
end
@kernel function tgeneric_kernel!(f!, v, args)
  i = @index(Global, Linear)
  f!(i, v, args...)
end

@kernel function tlindrft!(v, args)
  i = @index(Global, Linear)
  LinearTracking.linear_drift!(i, v, args...)
end
# collective effects
# each threads corresponds to many particles
# go through each element, each thread loops through each 
# particle and does stuff with it

# Call launch!
@inline runkernel!(f!::F, i::Nothing, v, work, args...) where {F} =launch!(f!, v, work, args...)

# Call kernel directly
@inline runkernel!(f!::F, i, v, work, args...) where {F} = f!(i, v, work, args...)


#=

for particle in particles
  for ele in ring

  end
end

for ele in ring
  # do a bunch pre pro
  for particle in particle

  end
end
 =#