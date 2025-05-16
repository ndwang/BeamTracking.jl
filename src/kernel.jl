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
  v::V, 
  args...; 
  groupsize::Union{Nothing,Integer}=nothing, #backend isa CPU ? floor(Int,REGISTER_SIZE/sizeof(eltype(v))) : 256 
  multithread_threshold::Integer=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
  use_KA::Bool=!(get_backend(v) isa CPU && isnothing(groupsize)),
  use_explicit_SIMD::Bool=false
) where {F<:Function,V}

  if use_KA && use_explicit_SIMD
    error("Cannot use both KernelAbstractions (KA) and explicit SIMD")
  end

  N_particle = size(v, 1)
  backend = get_backend(v)
  if !use_KA && backend isa GPU
    error("For GPU parallelized kernel launching, KernelAbstractions (KA) must be used")
  end

  if !use_KA
    if use_explicit_SIMD && V <: SIMD.FastContiguousArray && eltype(V) <: SIMD.ScalarTypes && VectorizationBase.pick_vector_width(eltype(V)) > 1 # do SIMD
      simd_lane_width = VectorizationBase.pick_vector_width(eltype(V))
      lane = VecRange{Int(simd_lane_width)}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          f!(lane+i, v, args...)
        end
      else
        for i in 1:simd_lane_width:N_SIMD
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          f!(lane+i, v, args...)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, args...)
      end
    else
      if N_particle >= multithread_threshold
        Threads.@threads for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          f!(i, v, args...)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          f!(i, v, args...)
        end
      end
    end
  else
    if isnothing(groupsize)
      kernel! = f!(backend)
    else
      kernel! = f!(backend, groupsize)
    end
    kernel!(v, args...; ndrange=N_particle)
    KernelAbstractions.synchronize(backend)
  end
  return v
end

# TODO: collective effects

# Helper functions for kernel execution
# When running kernel on a bunch, no index is provided, launch the kernel with automatic optimization
@inline runkernel!(f!::F, i::Nothing, v, args...; kwargs...) where {F} =launch!(f!, v, args...; kwargs...)
# When running kernel on a specific particle, execute the kernel directly for particle at that index
@inline runkernel!(f!::F, i, v, args...; kwargs...) where {F} = f!(i, v, args...)


macro makekernel(fcn)
  fcn.head == :function || error("@makekernel must wrap a function definition")
  body = esc(fcn.args[2])
  signature = fcn.args[1].args

  fcn_name = esc(signature[1])
  args = esc.(signature[2:end])
  i = esc(signature[2])
  v = esc(signature[3])
  work = esc(signature[4])

  const_args = map(signature[5:end]) do t
    if t isa Expr
      if t.head == :(::)
        :(@Const($(esc(t))))
      else
        error("Default values and keyword arguments are NOT supported by @Const in KernelAbstractions.jl")
      end
    else
      :(@Const($(esc(t))))
    end
  end
  stripped_args = map(signature[2:end]) do t
    if t isa Expr
      if t.args[1] isa Expr
        esc(t.args[1].args[1])
      else
        esc(t.args[1])
      end
    else
      esc(t)
    end
  end

  return quote
    @kernel function $(fcn_name)($v, $work, $(const_args...))
      $(stripped_args[1]) = @index(Global, Linear)
      $(fcn_name)($(stripped_args...))
    end
  
    @inline function $(fcn_name)($(args...))
        $(body)
    end
  end
end