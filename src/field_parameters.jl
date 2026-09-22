@inline _rebuild_functional_field(field_source::FunctionalField{F,P,N}, parameters) where {F,P,N} =
  FunctionalField(field_source.evaluator, parameters; normalized=N)

@inline batch_lower(field_source::FunctionalField) =
  _rebuild_functional_field(field_source, batch_lower(field_source.parameters))
@inline time_lower(field_source::FunctionalField) =
  _rebuild_functional_field(field_source, time_lower(field_source.parameters))
@inline num_lower(::Type{T}, field_source::FunctionalField) where {T<:Union{Float32,Float16}} =
  _rebuild_functional_field(field_source, num_lower(T, field_source.parameters))
@inline static_batchcheck(field_source::FunctionalField) =
  static_batchcheck(field_source.parameters)
@inline static_timecheck(field_source::FunctionalField) =
  static_timecheck(field_source.parameters)
@inline beval(field_source::FunctionalField, i) =
  _rebuild_functional_field(field_source, beval(field_source.parameters, i))
@inline teval(field_source::FunctionalField, t) =
  _rebuild_functional_field(field_source, teval(field_source.parameters, t))

@inline _rebuild_sum_field(field_source::SumField, field_sources) =
  SumField{typeof(field_sources)}(field_sources)

@inline batch_lower(field_source::SumField) =
  _rebuild_sum_field(field_source, batch_lower(field_source.field_sources))
@inline time_lower(field_source::SumField) =
  _rebuild_sum_field(field_source, time_lower(field_source.field_sources))
@inline num_lower(::Type{T}, field_source::SumField) where {T<:Union{Float32,Float16}} =
  _rebuild_sum_field(field_source, num_lower(T, field_source.field_sources))
@inline static_batchcheck(field_source::SumField) =
  static_batchcheck(field_source.field_sources)
@inline static_timecheck(field_source::SumField) =
  static_timecheck(field_source.field_sources)
@inline beval(field_source::SumField, i) =
  _rebuild_sum_field(field_source, beval(field_source.field_sources, i))
@inline teval(field_source::SumField, t) =
  _rebuild_sum_field(field_source, teval(field_source.field_sources, t))
