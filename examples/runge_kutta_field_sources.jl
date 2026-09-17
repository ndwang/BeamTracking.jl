using BeamTracking
using Beamlines
using Printf

# Run from an environment with BeamTracking and Beamlines installed:
#   julia --project=/path/to/environment examples/runge_kutta_field_sources.jl

const INITIAL_COORDS = [
   1.0e-3   1.0e-2  -2.0e-3   3.0e-3   0.0   0.0
  -1.5e-3  -4.0e-3   1.0e-3  -8.0e-3   0.0   1.0e-3
]

# StaticArrays is an implementation dependency of BeamTracking. Keeping this
# helper local means the demo environment only needs BeamTracking and Beamlines
# as direct dependencies.
svector(value) = BeamTracking.StaticArrays.SVector{1}(value)

"""A parameterized field evaluator suitable for `FunctionalField`."""
function uniform_magnetic_field(x, y, s, t, parameters)
  carrier = zero(x)
  return EMField(
    carrier, carrier, carrier,
    carrier + parameters.Bx,
    carrier + parameters.By,
    carrier + parameters.Bz,
  )
end

"""An example of a custom concrete callable field source."""
struct UniformBy{T}
  strength::T
end

function (field_source::UniformBy)(x, y, s, t)
  carrier = zero(x)
  return EMField(
    carrier, carrier, carrier,
    carrier, carrier + field_source.strength, carrier,
  )
end

function track_element(element, species, p_over_q_ref)
  line = Beamline(
    [element];
    species_ref=species,
    p_over_q_ref=p_over_q_ref,
  )
  bunch = Bunch(copy(INITIAL_COORDS); species=species, p_over_q_ref=p_over_q_ref)
  track!(bunch, line)
  return copy(bunch.coords.v)
end

function verify_case(label, actual, expected; atol=5e-13, rtol=5e-12)
  error_value = maximum(abs, actual - expected)
  passed = all(isapprox.(actual, expected; atol=atol, rtol=rtol))
  @printf("%-38s max |error| = %.3e  %s\n", label, error_value, passed ? "PASS" : "FAIL")
  passed || error("Runge-Kutta demo failed for: $label")
  return nothing
end

function main()
  species = Species("electron")
  p_over_q_ref = BeamTracking.pc_to_R(species, 1.0e9)
  length = 0.5
  n_steps = 20

  println("Runge-Kutta field-source demo")
  particle_count = size(INITIAL_COORDS, 1)
  println("Tracking $particle_count electrons at 1 GeV through $(length) m elements.\n")

  # 1. An explicit ZeroField under RK should reproduce exact drift tracking.
  exact_drift = track_element(Drift(L=length), species, p_over_q_ref)
  rk_drift = track_element(
    Drift(L=length, field_source=ZeroField(), tracking_method=RungeKutta(n_steps=n_steps)),
    species,
    p_over_q_ref,
  )
  verify_case("ZeroField vs exact drift", rk_drift, exact_drift)

  # 2. The element's normalized quadrupole and an explicit physical multipole
  #    describe the same magnetic field. Multipole order 2 is a quadrupole.
  kn1 = 0.2
  element_quadrupole = track_element(
    Quadrupole(L=length, Kn1=kn1, tracking_method=RungeKutta(n_steps=n_steps)),
    species,
    p_over_q_ref,
  )
  explicit_quadrupole = MultipoleField(
    svector(2),
    svector(kn1 * p_over_q_ref),
    svector(0.0),
  )
  source_quadrupole = track_element(
    Drift(
      L=length,
      field_source=explicit_quadrupole, tracking_method=RungeKutta(n_steps=n_steps),
    ),
    species,
    p_over_q_ref,
  )
  verify_case("MultipoleField vs element field", source_quadrupole, element_quadrupole)

  # 3. A FunctionalField can reproduce a uniform dipole field.
  uniform = FunctionalField(
    uniform_magnetic_field,
    (Bx=0.0, By=4.0e-3, Bz=0.0),
  )
  functional_result = track_element(
    Drift(L=length, field_source=uniform, tracking_method=RungeKutta(n_steps=n_steps)),
    species,
    p_over_q_ref,
  )
  dipole_result = track_element(
    Drift(
      L=length,
      field_source=MultipoleField(svector(1), svector(4.0e-3), svector(0.0)), tracking_method=RungeKutta(n_steps=n_steps),
    ),
    species,
    p_over_q_ref,
  )
  verify_case("FunctionalField vs uniform dipole", functional_result, dipole_result)

  # 4. additional_field composes an external field source with the element field.
  combined_element = track_element(
    Quadrupole(
      L=length,
      Kn1=kn1,
      additional_field=uniform, tracking_method=RungeKutta(n_steps=n_steps),
    ),
    species,
    p_over_q_ref,
  )
  combined_field_source = track_element(
    Drift(
      L=length,
      field_source=SumField(explicit_quadrupole, uniform), tracking_method=RungeKutta(n_steps=n_steps),
    ),
    species,
    p_over_q_ref,
  )
  verify_case("SumField vs additional_field", combined_field_source, combined_element)

  # 5. A concrete callable object can be used directly as a field source.
  custom_result = track_element(
    Drift(
      L=length,
      field_source=UniformBy(4.0e-3), tracking_method=RungeKutta(n_steps=n_steps),
    ),
    species,
    p_over_q_ref,
  )
  verify_case("Custom callable vs FunctionalField", custom_result, functional_result)

  println("\nAll field-source cases passed.")
  @printf(
    "Combined-field final coordinates (particle 1):\n  x=% .6e  px=% .6e  y=% .6e  py=% .6e  zeta=% .6e  delta=% .6e\n",
    combined_field_source[1, :]...,
  )
  return nothing
end

main()
