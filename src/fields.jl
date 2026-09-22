"""
    EMField{T}

Electric and magnetic field vectors at one particle location, in either physical
or normalized units as declared by the field source. In physical units, electric
field components are in V/m and magnetic field components are in tesla. In
normalized units, both vectors are divided by reference rigidity
`p_over_q_ref = p₀/q`, matching the four-potential convention.
`EMField` itself does not store a units flag or convert the supplied values.
"""
struct EMField{T}
  E::SVector{3,T}
  B::SVector{3,T}
end

function EMField(E::SVector{3,TE}, B::SVector{3,TB}) where {TE,TB}
  Ep, Bp = promote(E, B)
  return EMField{eltype(Ep)}(Ep, Bp)
end

function EMField(Ex, Ey, Ez, Bx, By, Bz)
  Ex, Ey, Ez, Bx, By, Bz = promote(Ex, Ey, Ez, Bx, By, Bz)
  return EMField(SVector(Ex, Ey, Ez), SVector(Bx, By, Bz))
end

@inline Base.:+(a::EMField, b::EMField) = EMField(a.E + b.E, a.B + b.B)

"""
    ZeroField()

A callable electromagnetic field source which returns zero fields.
Field sources are called as `field_source(x, y, s, t)`, with spatial coordinates in metres
and element-local particle time `t` in seconds (reference entrance at zero).
"""
struct ZeroField end

@inline function (::ZeroField)(x, y, s, t)
  v = zero(x)
  zero_vector = SVector(v, v, v)
  return EMField(zero_vector, zero_vector)
end

"""
    multipole_field(x, y, s, t, parameters)

Evaluate magnetic multipoles from `(orders, normal, skew)` coefficients.
Element tracking obtains these coefficients exclusively from `BMultipoleParams`
during unpacking. The returned field uses the coefficients' units.
"""
@inline function multipole_field(x, y, s, t, parameters)
  orders, normal, skew = parameters
  bx, by = normalized_field(orders, normal, skew, x, y, 0)
  zero_field = zero(bx)
  bz = vifelse(orders[1] == 0, normal[1], zero_field)
  return EMField(SVector(zero_field, zero_field, zero_field), SVector(bx, by, bz))
end

# Keep the evaluator and its payload explicit so all parameter preparation
# traverses the payload, including when the user's payload is nothing.
@inline call_field_function(x, y, s, t, parameters) =
  parameters[1](x, y, s, t, parameters[2])

"""
    FunctionalField(evaluator[, parameters]; normalized=false)

A callable field source backed by a concrete evaluator. With parameters, the
evaluator is called as `evaluator(x, y, s, t, parameters)`. Without parameters,
it is called as `evaluator(x, y, s, t)`. The evaluator must return `EMField`.
Spatial coordinates are in metres and `t` is in seconds, with zero at the
reference particle's element entrance. RK supplies each stage's particle time.
With `normalized=false` (the default), the evaluator must return E in V/m and B
in tesla. With `normalized=true`, it must return both E and B already divided by
reference rigidity `p_over_q_ref = p₀/q`. The flag declares the evaluator's output
units; direct calls preserve those values. Tracking converts physical fields to
normalized units and uses normalized fields directly.
"""
struct FunctionalField{F,P,N}
  evaluator::F
  parameters::P
end

FunctionalField(evaluator, parameters=nothing; normalized::Bool=false) =
  FunctionalField{typeof(evaluator),typeof(parameters),normalized}(evaluator, parameters)

@inline function (field_source::FunctionalField{F,Nothing,N})(x, y, s, t) where {F,N}
  return field_source.evaluator(x, y, s, t)
end

@inline function (field_source::FunctionalField)(x, y, s, t)
  return field_source.evaluator(x, y, s, t, field_source.parameters)
end

"""
    SumField(field_sources...)
    SumField(field_sources::Tuple)

A callable, statically dispatched sum of electromagnetic field sources.
Nested sums are flattened and `ZeroField` members are removed at construction.
Direct evaluation requires all components to use the same units. Tracking
normalizes each component before addition and also supports mixed-unit sums.
"""
struct SumField{S<:Tuple}
  field_sources::S

  SumField{S}(field_sources::S) where {S<:Tuple} = new{S}(field_sources)
end

@inline _flatten_field_source(::ZeroField) = ()
@inline _flatten_field_source(field_source::SumField) = field_source.field_sources
@inline _flatten_field_source(field_source) = (field_source,)

@inline _flatten_field_sources(::Tuple{}) = ()
@inline function _flatten_field_sources(field_sources::Tuple)
  return (_flatten_field_source(first(field_sources))..., _flatten_field_sources(Base.tail(field_sources))...)
end

function SumField(field_sources::Tuple)
  flattened = _flatten_field_sources(field_sources)
  isempty(flattened) && return ZeroField()
  length(flattened) == 1 && return first(flattened)
  return SumField{typeof(flattened)}(flattened)
end

SumField(field_sources...) = SumField(field_sources)

@generated function _evaluate_field_sum(field_sources::S, x, y, s, t) where {S<:Tuple}
  N = length(S.parameters)
  N > 0 || return :(ZeroField()(x, y, s, t))
  expression = :(Base.getfield(field_sources, 1)(x, y, s, t))
  for i in 2:N
    expression = :($expression + Base.getfield(field_sources, $i)(x, y, s, t))
  end
  return expression
end

# Unit traits are compile-time constants, like the implicit integrator's Val flag.
@inline field_normalized(field_source) = Val(false)
@inline field_normalized(::FunctionalField{F,P,N}) where {F,P,N} = Val(N)
@inline function field_normalized(field_source::SumField)
  units = field_normalized(first(field_source.field_sources))
  all(s -> field_normalized(s) == units, field_source.field_sources) ||
    throw(ArgumentError("mixed-unit SumField requires reference rigidity; evaluate it through tracking"))
  return units
end

"""
    normalized_field_at(field_source, x, y, s, t, inv_rigidity)

Evaluate fields with both E and B divided by reference rigidity, as in
`implicit_fields`. Custom field sources default to physical units; wrap a normalized
custom evaluator in `FunctionalField(...; normalized=true)`.
"""
@inline normalized_field_at(field_source, x, y, s, t, inv_rigidity) =
  normalized_field_at(field_source, x, y, s, t, inv_rigidity, field_normalized(field_source))

@inline function normalized_field_at(field_source, x, y, s, t, inv_rigidity, ::Val{normalized}) where {normalized}
  field = field_source(x, y, s, t)
  if normalized
    return field
  else
    return EMField(field.E * inv_rigidity, field.B * inv_rigidity)
  end
end

@inline function normalized_field_at(field_source::SumField, x, y, s, t, inv_rigidity)
  fields = map(field_source.field_sources) do f
    @inline
    normalized_field_at(f, x, y, s, t, inv_rigidity)
  end
  return +(fields...)
end

@inline function (field_source::SumField)(x, y, s, t)
  field_normalized(field_source) # Reject adding physical and normalized values directly.
  return _evaluate_field_sum(field_source.field_sources, x, y, s, t)
end

include("field_parameters.jl")

Adapt.@adapt_structure EMField
@inline function Adapt.adapt_structure(to, field_source::FunctionalField{F,P,N}) where {F,P,N}
  return FunctionalField(Adapt.adapt(to, field_source.evaluator),
                         Adapt.adapt(to, field_source.parameters); normalized=N)
end
Adapt.@adapt_structure SumField
