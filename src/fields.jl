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
    zero_field(x, y, s, t, parameters=nothing)

Return zero electric and magnetic fields with the coordinate's numeric type.
"""
@inline function zero_field(x, y, s, t, parameters=nothing)
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

"""
    normalized_field_at(field_function, parameters, normalized, x, y, s, t, inv_rigidity)

Evaluate `field_function(x, y, s, t, parameters)` and return both E and B
normalized by reference rigidity. `normalized` is `Val(true)` for an evaluator
that already returns normalized fields, and `Val(false)` for physical fields.

Parallel tuples of functions, parameters, and flags add contributions after
normalizing each independently. Empty tuples produce zero field.
"""
@inline function normalized_field_at(field_function, parameters, ::Val{normalized},
                                     x, y, s, t, inv_rigidity) where {normalized}
  field = field_function(x, y, s, t, parameters)
  if normalized
    return field
  else
    return EMField(field.E * inv_rigidity, field.B * inv_rigidity)
  end
end

@inline function normalized_field_at(field_functions::Tuple, parameters::Tuple, normalized::Tuple,
                                     x, y, s, t, inv_rigidity)
  return _normalized_field_sum(field_functions, parameters, normalized, x, y, s, t, inv_rigidity)
end

@inline _normalized_field_sum(::Tuple{}, ::Tuple{}, ::Tuple{}, x, y, s, t, inv_rigidity) =
  zero_field(x, y, s, t)

@inline function _normalized_field_sum(field_functions::Tuple, parameters::Tuple, normalized::Tuple,
                                       x, y, s, t, inv_rigidity)
  field = normalized_field_at(first(field_functions), first(parameters), first(normalized),
                              x, y, s, t, inv_rigidity)
  # Avoid adding a zero field to the last contribution, especially for TPSA.
  length(field_functions) == 1 && return field
  return field + _normalized_field_sum(Base.tail(field_functions), Base.tail(parameters),
                                       Base.tail(normalized), x, y, s, t, inv_rigidity)
end

Adapt.@adapt_structure EMField
