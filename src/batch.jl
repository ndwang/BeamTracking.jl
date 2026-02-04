#= 

BatchParam is very similar to TimeDependentParam, but instead 
of storing a function, it stores an arbitrary array of parameters.

Given a batch = [k1, k2, k3], batch parameters are seen by the particles as

Particle 1: k1
Particle 2: k2
Particle 3: k3
Particle 4: k1
Particle 5: k2
Particle 6: k3

etc.

Just like time, there are two types for batches - one type unstable 
generic wrapper for an AbstractArray/Number which is manipulated at the 
highest level, and a lowered representation where the numbers are made
numbers and arrays are kept but the length of the array is stored in 
the type. The reason for this is so that scalars are not unnecessarily
represented as large arrays, to save memory both outside and inside the 
kernel. The goal is that for large, batched simulations the type instability
of the UNPACKING step (simulation step is always type stable) is outweighted 
by the simulation step.

The lowered type does NOT have any arithmetic operations defined on it, as 
it should be untouched after lowering and in the kernel evaluated for each 
particle.

=#
struct BatchParam
  batch::Union{AbstractArray,Number}
  #=
    Currently, batches of TimeDependentParam are not supported
    It is essentially impossible on the GPU because accessing an 
    array of functions is not type stable (every particle has its 
    own function). The alternative - a time function which returns a
    big array - is also not possible for CPU-SIMD and would run into 
    memory problems on both CPU and GPU depending on how large that 
    static array is.
    
    On the CPU, the first option may be doable using FunctionWrappers,
    but that would require special handling because current Time
    uses bona-fide Julia functions for GPU compatibility.

    If one would like to do this right now, they should just start up 
    separate processes where each process has its own lattice with its 
    own TimeDependentParams.
  =#
end

struct _LoweredBatchParam{N,V<:AbstractArray}
  batch::V
  _LoweredBatchParam(batch::AbstractArray) = new{length(batch),typeof(batch)}(batch)
  _LoweredBatchParam{N}(batch::AbstractArray) where {N} = new{N,typeof(batch)}(batch)
end

# Necessary for GPU compatibility if batch is a GPU array
function Adapt.adapt_structure(to, lbp::_LoweredBatchParam{N}) where {N}
    batch = Adapt.adapt_structure(to, lbp.batch)
    return _LoweredBatchParam{N,typeof(batch)}(batch)
end

# BatchParam will act like a number
# Conversion of types to BatchParam
BatchParam(a::BatchParam) = a

# Make these apply via convert
Base.convert(::Type{BatchParam}, a::Number) = BatchParam(a) # Scalar BatchParam
Base.convert(::Type{BatchParam}, a::BatchParam) = a

Base.zero(b::BatchParam) = BatchParam(zero(first(b.batch)))
Base.one(b::BatchParam)  = BatchParam(one(first(b.batch))) 

# Now define the math operations:
for op in (:+,:-,:*,:/,:^)
  @eval begin
    Base.$op(ba::BatchParam, b::Number)   = (let b = b; return BatchParam(map(x->$op(x, b), ba.batch)); end)
    Base.$op(a::Number,   bb::BatchParam) = (let a = a; return BatchParam(map(x->$op(a, x), bb.batch)); end)
    function Base.$op(ba::BatchParam, bb::BatchParam)
      if !(ba.batch isa Number) && !(bb.batch isa Number) && length(ba.batch) != length(bb.batch)
        error("Cannot perform operation $($op) with two non-scalar BatchParams of differing 
               lengths (received lengths $(length(ba.batch)) and $(length(bb.batch))).")
      end
      return BatchParam(map((x,y)->$op(x, y), ba.batch, bb.batch))
    end
  end
end

function Base.literal_pow(::typeof(^), ba::BatchParam, ::Val{N}) where {N}
  return BatchParam(map(x->Base.literal_pow(^, x, Val{N}()), ba.batch))
end

for t = (:+, :-, :sqrt, :exp, :log, :sin, :cos, :tan, :cot, :sinh, :cosh, :tanh, :inv,
  :coth, :asin, :acos, :atan, :acot, :asinh, :acosh, :atanh, :acoth, :sinc, :csc, :float,
  :csch, :acsc, :acsch, :sec, :sech, :asec, :asech, :conj, :log10, :isnan, :sign, :abs)
  @eval begin
    Base.$t(b::BatchParam) = BatchParam(map(x->($t)(x), b.batch))
  end
end

atan2(b1::BatchParam, b2::BatchParam) = BatchParam(map((x,y)->atan2(x,y), b1.batch, b2.batch))

for t = (:unit, :sincu, :sinhc, :sinhcu, :asinc, :asincu, :asinhc, :asinhcu, :erf, 
         :erfc, :erfcx, :erfi, :wf, :rect)
  @eval begin
    GTPSA.$t(b::BatchParam) = BatchParam(map(x->($t)(x), b.batch))
  end
end


Base.promote_rule(::Type{BatchParam}, ::Type{U}) where {U<:Number} = BatchParam
Base.broadcastable(o::BatchParam) = Ref(o)

Base.isapprox(b::BatchParam, n::Number; kwargs...) = all(x->isapprox(x, n, kwargs...), b.batch)
Base.isapprox(n::Number, b::BatchParam; kwargs...) = all(x->isapprox(n, x, kwargs...), b.batch)
Base.:(==)(b::BatchParam, n::Number) = all(x->x == n, b.batch)
Base.:(==)(n::Number, b::BatchParam) = all(x->n == x, b.batch)
Base.isinf(b::BatchParam) = all(x->isinf(x), b.batch)

#=
@inline teval(f::TimeFunction, t) = f(t)
@inline teval(f, t) = f

# === THIS BLOCK WAS WRITTEN BY CLAUDE ===
# Generated function for arbitrary-length tuples
@generated function teval(f::T, t) where {T<:Tuple}
    N = length(T.parameters)
    if N == 0
        return :(())
    end
    # Use getfield with literal integer arguments
    exprs = [:(teval(Base.getfield(f, $i), t)) for i in 1:N]
    return :(tuple($(exprs...)))
end
# === END CLAUDE ===
=#

# Batch lowering should convert types to _LoweredBatchParam
function batch_lower(b::BatchParam)
  if b.batch isa AbstractArray
    return _LoweredBatchParam(b.batch)
  else
    return b.batch
  end
end

batch_lower(tp::TimeDependentParam) = tp.f
batch_lower(tp) = tp
# We can use map on the CPU, but not the GPU. This step of time_lower-ing is on 
# the CPU and we are already type unstable here anyways, so we should do this.
time_lower(tp::T) where {T<:Tuple} = map(ti->time_lower(ti), tp)

# Arrays MUST be converted into tuples, for SIMD
time_lower(tp::SArray{N,TimeDependentParam}) where {N} = time_lower(Tuple(tp))

static_timecheck(tp) = false
static_timecheck(::TimeFunction) = true
@unroll function static_timecheck(t::Tuple)
  @unroll for ti in t
    if static_timecheck(ti)
      return true
    end
  end
  return false
end
#static_timecheck(tp::T) where {T<:Tuple} = Val{any(t->static_timecheck(t) isa Val{true}, tp)}()



batch_lower(b) = b
batch_lower(b::Tuple) = map(bi->batch_lower(bi), b)
function batch_lower(ba::SArray{N,BatchParam}) where {N}
  f = Tuple(map(bi->bi.batch, ba))
  return TimeFunction(t->SArray{N}(map(fi->fi(t), f)))
end
time_lower(tp::SArray{N,Any}) where {N} = time_lower(TimeDependentParam.(tp))

static_timecheck(tp) = false
static_timecheck(::TimeFunction) = true
@unroll function static_timecheck(t::Tuple)
  @unroll for ti in t
    if static_timecheck(ti)
      return true
    end
  end
  return false
end

"""
    lane2vec(lane::SIMD.VecRange{N}) 
    
Given a SIMD.VecRange, will return an equivalent SIMD.Vec that
can be used in arithmetic operations for mapping integer indices
of particles to a given element in a batch.
"""
function lane2vec(lane::SIMD.VecRange{N}) where {N}
  # Try to match with vector register size, but 
  # only up to UInt32 -> ~4.3 billion particles, 
  # probably max on CPU...
  if Int(pick_vector_width(UInt32)) == N
    return SIMD.Vec{N,UInt32}(ntuple(i->lane.i+i-1, Val{N}()))
  else
    return SIMD.Vec{N,UInt64}(ntuple(i->lane.i+i-1, Val{N}()))
  end
end

@inline teval(f::TimeFunction, t) = f(t)
@inline teval(f, t) = f
@inline teval(f::Tuple, t) = map(ti->teval(ti, t), f)
