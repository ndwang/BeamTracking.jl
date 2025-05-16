# Get the register size for SIMD operations from VectorizationBase
const REGISTER_SIZE = VectorizationBase.register_size()

"""
    launch!(f!::F, v::A, work, args...; simd_lane_width, multithread_threshold)

Launch a kernel function on particle coordinates with automatic optimization.

# Arguments
- `f!`: Kernel function to execute. The kernel function `f!` must be of the form `f!(i, v, work, args...)`
- `v`: Input/output matrix of particle coordinates (always in SoA format)
- `work`: Vector of temporary vectors for kernel execution
- `args...`: Additional arguments for the kernel function

# Keyword Arguments
- `simd_lane_width`: Number of SIMD lanes to use. Default is 0 (autovectorize)
- `multithread_threshold`: Particle count threshold for multithreading (default: 1750 * nthreads)
"""
@inline function launch!(
  f!::F, 
  v::A,
  work, 
  args...; 
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {F<:Function,A}
  N_particle = size(v, 1)
  
  # SIMD path: Use vectorized operations when possible
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0
    # Create a vector range for SIMD operations
    lane = VecRange{simd_lane_width}(0)
    # Calculate remainder for non-SIMD processing
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    
    # Multithreaded SIMD path
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    # Single-threaded SIMD path
    else
      for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    end
    
    # Process remaining particles that don't fit in SIMD lanes
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, work, args...)
    end
  # Non-SIMD path: Use standard loops with potential multithreading
  else
    if N_particle >= multithread_threshold
      # Multithreaded standard path
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    else
      # Single-threaded standard path with @simd hint
      @simd for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    end
  end
  return v
end

# TODO: collective effects

# Helper functions for kernel execution
# When running kernel on a bunch, no index is provided, launch the kernel with automatic optimization
@inline runkernel!(f!::F, i::Nothing, v, work, args...) where {F} = launch!(f!, v, work, args...)
# When running kernel on a specific particle, execute the kernel directly for particle at that index
@inline runkernel!(f!::F, i, v, work, args...) where {F} = f!(i, v, work, args...)
