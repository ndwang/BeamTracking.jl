# The TimeDependentParam is quite different from a 
# deferred expression as in Beamlines. For this, there
# are actually no closures used - the TimeDependentParam
# is constructed once the deferred expression is deval-ed.
# Therefore not nesting of the TimeDependentParams to 
# capture the "state" is necessary - they will all be just 
# Functions of regular numbers

struct TimeDependentParam
  f::Function
  TimeDependentParam(f::Function=(t)->t) = new(f)
end

# Convenience ctor
Time() = TimeDependentParam()

# Calling TimeDependentParam
(d::TimeDependentParam)(t) = d.f(t)

# Conversion of types to TimeDependentParam
TimeDependentParam(a::Number) = TimeDependentParam((t)->a) 
TimeDependentParam(a::TimeDependentParam) = a

# Make these apply via convert
Base.convert(::Type{T}, a) where {D<:TimeDependentParam} = D(a)

#teval(tp::TimeDependentParam, t) = tp(t)
#teval(tp, t) = tp

# Now define the math operations:
for op in (:+,:-,:*,:/,:^)
  @eval begin
    Base.$op(da::TimeDependentParam, b)   = (let fa = da.f, b = b; return TimeDependentParam((t)-> $op(fa(t), b)); end)
    Base.$op(a,   db::TimeDependentParam) = (let fb = db.f, a = a; return TimeDependentParam((t)-> $op(a, fb(t))); end)
    function Base.$op(da::TimeDependentParam, db::TimeDependentParam)
      let fa = da.f, fb = db.f
        return TimeDependentParam((t)-> $op(fa(t), fb(t)))
      end
    end
  end
end

for t = (:sqrt, :exp, :log, :sin, :cos, :tan, :cot, :sinh, :cosh, :tanh, :inv,
  :coth, :asin, :acos, :atan, :acot, :asinh, :acosh, :atanh, :acoth, :sinc, :csc, 
  :csch, :acsc, :acsch, :sec, :sech, :asec, :asech, :conj, :log10, :isnan)
  @eval begin
    Base.$t(d::TimeDependentParam) = (let f = d.f; return TimeDependentParam((t)-> ($t)(f(t))); end)
  end
end

for t = (:unit, :sincu, :sinhc, :sinhcu, :asinc, :asincu, :asinhc, :asinhcu, :erf, 
         :erfc, :erfcx, :erfi, :wf, :rect)
  @eval begin
    GTPSA.$t(d::TimeDependentParam) = (let f = d.f; return TimeDependentParam((t)-> ($t)(f(t))); end)
  end
end

Base.promote_rule(::Type{TimeDependentParam}, ::Type{U}) where {U<:Number} = TimeDependentParam
Base.broadcastable(o::TimeDependentParam) = Ref(o)