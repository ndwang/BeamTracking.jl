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
      if ba.batch isa Number
        if bb.batch isa Number
          return BatchParam($op(ba.batch, bb.batch))
        else
          let a = ba.batch
            # WE HAVE TO WRITE THIS BY HAND FOR EACH OP!!!
            if a ≈ 0 # So that e.g. big arrays aren't made for skew strengths when only normal is batch
              return BatchParam(0f0)
            else
              return BatchParam(map((bbi)->$op(a, bbi), bb.batch))
            end
          end
        end
      elseif bb.batch isa Number
        let b = bb.batch
          if b ≈ 0
            return BatchParam(0f0)
          else
            return BatchParam(map((bai)->$op(bai, b), ba.batch))
          end
        end
      elseif length(ba.batch) == length(bb.batch)
        return BatchParam(map((bai,bbi)->$op(bai, bbi), ba.batch, bb.batch))
      else
        error("Cannot perform operation $($op) with two non-scalar BatchParams of differing 
               lengths (received lengths $(length(ba.batch)) and $(length(bb.batch))).")
      end
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
Base.promote_rule(::Type{BatchParam}, ::Type{TimeDependentParam}) = error("Unable to combine BatchParams with TimeDependentParams")
Base.promote_rule(::Type{TimeDependentParam}, ::Type{BatchParam}) = error("Unable to combine BatchParams with TimeDependentParams")
Base.broadcastable(o::BatchParam) = Ref(o)

Base.isapprox(b::BatchParam, n::Number; kwargs...) = all(x->isapprox(x, n, kwargs...), b.batch)
Base.isapprox(n::Number, b::BatchParam; kwargs...) = all(x->isapprox(n, x, kwargs...), b.batch)
Base.:(==)(b::BatchParam, n::Number) = all(x->x == n, b.batch)
Base.:(==)(n::Number, b::BatchParam) = all(x->n == x, b.batch)
Base.isinf(b::BatchParam) = all(x->isinf(x), b.batch)

# Batch lowering should convert types to _LoweredBatchParam
function batch_lower(b::BatchParam)
  if b.batch isa AbstractArray
    return _LoweredBatchParam(b.batch) # Only arrays are lowered to batchparams
  else
    return b.batch
  end
end

batch_lower(bp) = bp
# We can use map on the CPU, but not the GPU. This step of batch_lower-ing is on 
# the CPU and we are already type unstable here anyways, so we should do this.
batch_lower(bp::T) where {T<:Tuple} = map(bi->batch_lower(bi), bp)

# Arrays MUST be converted into tuples, for SIMD
batch_lower(bp::SArray{N,BatchParam}) where {N} = batch_lower(Tuple(bp))

static_batchcheck(bp) = false
static_batchcheck(::BatchParam) = true
@unroll function static_batchcheck(t::Tuple)
  @unroll for ti in t
    if static_batchcheck(ti)
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

lane2vec(i) = i

@inline beval(b::_LoweredBatchParam, i) = b.batch[lane2vec(i)]
@inline beval(b, i) = b

# === THIS BLOCK WAS WRITTEN BY CLAUDE ===
# Generated function for arbitrary-length tuples
@generated function beval(f::T, t) where {T<:Tuple}
    N = length(T.parameters)
    if N == 0
        return :(())
    end
    # Use getfield with literal integer arguments
    exprs = [:(beval(Base.getfield(f, $i), t)) for i in 1:N]
    return :(tuple($(exprs...)))
end
# === END CLAUDE ===