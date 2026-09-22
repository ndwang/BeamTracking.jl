function test_uniform_field(x, y, s, t, parameters)
  carrier = zero(x)
  return EMField(
    carrier + parameters.Ex, carrier + parameters.Ey, carrier + parameters.Ez,
    carrier + parameters.Bx, carrier + parameters.By, carrier + parameters.Bz,
  )
end

function test_parameter_free_field(x, y, s, t, parameters)
  @assert isnothing(parameters)
  carrier = zero(x)
  return EMField(carrier, carrier, carrier, carrier, carrier, carrier + 1)
end

struct FieldTestAdaptor end
function BeamTracking.Adapt.adapt_storage(::FieldTestAdaptor, values::Vector)
  return SVector{length(values)}(values)
end

struct FieldTestParameters{T}
  strength::T
end

@testset "Field functions" begin
  @testset "EMField" begin
    field = EMField(SA[1.0, 2.0, 3.0], SA[4.0, 5.0, 6.0])
    @test field.E == SA[1.0, 2.0, 3.0]
    @test field.B == SA[4.0, 5.0, 6.0]
    @test EMField(1, 2, 3, 4, 5, 6) == EMField(SA[1, 2, 3], SA[4, 5, 6])
    @test field + field == EMField(SA[2.0, 4.0, 6.0], SA[8.0, 10.0, 12.0])
  end

  @testset "Zero field" begin
    field = @inferred BeamTracking.zero_field(1.0, 2.0, 3.0, 4.0, nothing)
    @test field == EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 0.0])
    dual = ForwardDiff.Dual(1.0, 1.0)
    dual_field = @inferred BeamTracking.zero_field(dual, dual, dual, 0.0, nothing)
    @test eltype(dual_field.E) === typeof(dual)
    @test eltype(dual_field.B) === typeof(dual)
    @test @inferred(BeamTracking.normalized_field_at((), (), (), 1.0, 2.0, 3.0, 4.0, 0.5)) == field
  end

  @testset "Multipole field evaluator" begin
    solenoid = (SA[0], SA[1.5], SA[0.0])
    dipole = (SA[1], SA[2.0], SA[3.0])
    quadrupole = (SA[2], SA[4.0], SA[5.0])
    @test @inferred(BeamTracking.multipole_field(0.2, 0.3, 0.0, 0.0, solenoid)) ==
      EMField(SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.5])
    @test @inferred(BeamTracking.multipole_field(0.2, 0.3, 0.0, 0.0, dipole)) ==
      EMField(SA[0.0, 0.0, 0.0], SA[3.0, 2.0, 0.0])
    @test @inferred(BeamTracking.multipole_field(0.2, 0.3, 0.0, 0.0, quadrupole)) ==
      EMField(SA[0.0, 0.0, 0.0], SA[2.2, -0.7, 0.0])
    simd_x = SIMD.Vec{2,Float64}((0.2, 0.4))
    simd_y = SIMD.Vec{2,Float64}((0.3, 0.1))
    simd_field = @inferred BeamTracking.multipole_field(simd_x, simd_y, zero(simd_x), zero(simd_x), quadrupole)
    @test all(isapprox.(Tuple(simd_field.B[1]), (2.2, 2.4)))
    @test all(isapprox.(Tuple(simd_field.B[2]), (-0.7, 1.1)))
    descriptor = Descriptor(2, 1)
    tpsa_x, tpsa_y = @vars(descriptor)
    tpsa_field = @inferred BeamTracking.multipole_field(tpsa_x, tpsa_y, zero(tpsa_x), zero(tpsa_x), quadrupole)
    @test GTPSA.jacobian(collect(tpsa_field.B[1:2])) ≈ [5.0 4.0; 4.0 -5.0]
  end

  @testset "Parameter preparation" begin
    context = Beamlines.Context(strength=0.25)
    parameters = (strength=Beamlines.DefExpr{Float64}(c -> 2 * c.strength),
                  constants=(label="map", count=2))
    field_function = (x, y, s, t, p) -> begin
      v = zero(x)
      EMField(v, v, v, v, v + p.strength, v)
    end
    group = FieldFunctionParams(field_function=field_function, field_function_params=parameters)
    prepared = @inferred Beamlines.deval(group, context)
    @test prepared.field_function === field_function
    @test prepared.field_function_params.strength == 0.5
    @test prepared.field_function_params.constants === parameters.constants
    @test group.field_function_params.strength isa Beamlines.DefExpr
    context.strength = 0.75
    @test Beamlines.deval(group, context).field_function_params.strength == 1.5
    opaque = FieldTestParameters(parameters.strength)
    opaque_group = FieldFunctionParams(field_function=field_function, field_function_params=opaque)
    @test Beamlines.deval(opaque_group, context).field_function_params === opaque
    dual_group = FieldFunctionParams(field_function=field_function,
      field_function_params=(strength=ForwardDiff.Dual(0.5, 1.0), constants=parameters.constants))
    scalar_group = @inferred Beamlines.scalarize(dual_group)
    @test scalar_group.field_function === field_function
    @test scalar_group.field_function_params.strength === 0.5

    dynamic = (strength=BatchParam([0.5, 1.5]), constants=parameters.constants)
    lowered = BeamTracking.batch_lower(dynamic)
    @test @inferred(BeamTracking.static_batchcheck(lowered))
    selected = @inferred BeamTracking.beval(lowered, 2, 1)
    @test field_function(0.0, 0.0, 0.0, 0.0, selected).B[2] == 1.5
    @test_opt BeamTracking.beval(lowered, 2, 1)
    @test @ballocated(BeamTracking.beval($lowered, 2, 1)) == 0
    timed = (strength=2 * Time(), constants=parameters.constants)
    lowered_time = BeamTracking.time_lower(timed)
    @test @inferred(BeamTracking.static_timecheck(lowered_time))
    evaluated = @inferred BeamTracking.teval(lowered_time, 0.25)
    @test field_function(0.0, 0.0, 0.0, 0.0, evaluated).B[2] == 0.5
    @test @ballocated(BeamTracking.teval($lowered_time, 0.25)) == 0
    @test @inferred(test_parameter_free_field(0.0, 0.0, 0.0, 0.0, nothing)).B == SA[0.0, 0.0, 1.0]
  end

  @testset "Named-tuple batch offsets" begin
    values = [10.0, 20.0, 30.0]
    parameters = (field=(By=BatchParam(values),), scale=2.0)
    lowered = BeamTracking.batch_lower(parameters)
    for batch_start in (2, 4), i in 1:4
      selected = @inferred BeamTracking.beval(lowered, i, batch_start)
      @test selected.field.By == values[mod1(i + batch_start - 1, length(values))]
      @test selected.scale == 2.0
    end
    if !(VERSION < v"1.11" && Sys.ARCH == :x86_64)
      selected = @inferred BeamTracking.beval(lowered, SIMD.VecRange{4}(1), 2)
      @test Tuple(selected.field.By) == (20.0, 30.0, 10.0, 20.0)
      @test selected.scale == 2.0
    end
  end

  @testset "Multipole parameter preparation" begin
    parameters = (SA[1], SA[BatchParam([2.0, 3.0])], SA[BatchParam(0.0)])
    lowered = BeamTracking.batch_lower(parameters)
    @test BeamTracking.static_batchcheck(lowered)
    for i in 1:2
      selected = @inferred BeamTracking.beval(lowered, i, 1)
      @test @inferred(BeamTracking.multipole_field(0.0, 0.0, 0.0, 0.0, selected)).B == SA[0.0, i + 1.0, 0.0]
    end
    timed = (SA[1], SA[4.0 * Time()], SA[TimeDependentParam(0.0)])
    evaluated = @inferred BeamTracking.teval(BeamTracking.time_lower(timed), 0.5)
    @test @inferred(BeamTracking.multipole_field(0.0, 0.0, 0.0, 0.0, evaluated)).B == SA[0.0, 2.0, 0.0]
  end

  @testset "Additive evaluation and allocations" begin
    functions = (BeamTracking.multipole_field, test_uniform_field)
    parameters = ((SA[2], SA[4.0], SA[5.0]),
      (Ex=0.0, Ey=0.0, Ez=0.0, Bx=1.0, By=0.0, Bz=3.0))
    flags = (Val(true), Val(false))
    result = @inferred BeamTracking.normalized_field_at(functions, parameters, flags, 0.2, 0.3, 0.0, 0.0, 0.5)
    @test result.B ≈ SA[2.7, -0.7, 1.5]
    @test_opt BeamTracking.normalized_field_at(functions, parameters, flags, 0.2, 0.3, 0.0, 0.0, 0.5)
    @test @ballocated(BeamTracking.normalized_field_at($functions, $parameters, $flags, 0.2, 0.3, 0.0, 0.0, 0.5)) == 0
  end

  @testset "Numeric lowering" begin
    functions = (BeamTracking.multipole_field, test_uniform_field)
    multipole = (SA[1, 2], SA[0.01, 0.03], SA[0.0, 0.02])
    functional = (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0, By=0.1, Bz=0.0,
      nested=(values=(0.25, SA[0.5, 0.75]), order=2, label="field"))
    parameters = (multipole, functional)
    for T in (Float32, Float16)
      lowered = @inferred BeamTracking.num_lower(T, parameters)
      @test lowered[1][1] === multipole[1]
      @test lowered[2].nested.values === (T(0.25), SVector{2,T}(0.5, 0.75))
      @test lowered[2].nested.order === 2
      @test lowered[2].nested.label === "field"
      @test @inferred(BeamTracking.normalized_field_at(functions, lowered, (Val(true), Val(true)), zero(T), zero(T), zero(T), zero(T), one(T))) isa EMField{T}
    end
    @test BeamTracking.num_lower(Float64, parameters) === parameters
    opaque = FieldTestParameters(0.1)
    @test BeamTracking.num_lower(Float32, opaque) === opaque
    batch = (SA[1], SA[BatchParam([0.1, 0.2])], SA[BatchParam(0.0)])
    timed = (Ex=0.0, Ey=0.0, Ez=0.0, Bx=0.0,
      By=TimeDependentParam(t -> Float32(t), false), Bz=0.0)
    call = BeamTracking.make_kernel_call(identity, ((batch, timed),))
    chain = BeamTracking.KernelChain(Val{1}(), BeamTracking.RefState{Float32}(;
      t_enter=0f0, beta_gamma_enter=1f0))
    prepared = BeamTracking.push(chain, call).chain[1].args[1]
    for i in 1:2
      evaluated = BeamTracking.teval(BeamTracking.beval(prepared, i, 1), 0.25f0)
      field = @inferred BeamTracking.normalized_field_at(functions, evaluated, (Val(true), Val(true)), 0f0, 0f0, 0f0, 0f0, 1f0)
      @test field isa EMField{Float32}
      @test field.B[2] ≈ Float32(i * 0.1) + 0.25f0
    end
    @test batch[2][1].batch == [0.1, 0.2]
    bad_call = BeamTracking.make_kernel_call(identity, ((By=Time(),),))
    @test_throws ErrorException BeamTracking.push(chain, bad_call)
  end

  @testset "Adaptation" begin
    multipole = (SA[0, 2], SA[0.01, 0.03], SA[0.0, 0.02])
    @test BeamTracking.Adapt.adapt(FieldTestAdaptor(), multipole) == multipole
    @test @ballocated(BeamTracking.Adapt.adapt(FieldTestAdaptor(), $multipole)) == 0
    parameters = (multipole, (field_map=[1.0, 2.0, 3.0],))
    adapted = BeamTracking.Adapt.adapt(FieldTestAdaptor(), parameters)
    @test adapted[1] === multipole
    @test adapted[2].field_map == SA[1.0, 2.0, 3.0]
  end
end

@testset "Field unit conventions" begin
  for T in (Float32, Float64), R in (T(-3), T(2))
    args = (T(0.02), T(-0.01), zero(T), zero(T))
    field_function = (x,y,s,t,p) -> EMField(p...)
    physical_parameters = T.((1,2,3,4,5,6))
    expected = field_function(args..., physical_parameters)
    converted = @inferred BeamTracking.normalized_field_at(field_function, physical_parameters, Val(false), args..., inv(R))
    direct = @inferred BeamTracking.normalized_field_at(field_function, physical_parameters ./ R, Val(true), args..., inv(R))
    @test converted.E ≈ expected.E / R
    @test converted.B ≈ expected.B / R
    @test converted.E ≈ direct.E
    @test converted.B ≈ direct.B
    @test field_function(args..., physical_parameters).E == expected.E
    @test @inferred(BeamTracking.normalized_field_at(test_parameter_free_field, nothing, Val(true), args..., inv(R))) == test_parameter_free_field(args..., nothing)
    multipole = (SA[1,2], T.(SA[0.1,0.2]), T.(SA[0,0]))
    result = @inferred BeamTracking.normalized_field_at((BeamTracking.multipole_field, field_function), (multipole, physical_parameters), (Val(true), Val(false)), args..., inv(R))
    @test result.E ≈ converted.E
    @test result.B ≈ BeamTracking.multipole_field(args..., multipole).B + converted.B
  end
end
