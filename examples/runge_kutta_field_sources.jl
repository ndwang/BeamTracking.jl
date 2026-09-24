using BeamTracking
using Beamlines
using Printf

# Run from an environment with BeamTracking and Beamlines installed:
#   julia --project=/path/to/environment examples/runge_kutta_field_sources.jl

const INITIAL_COORDS = [
   1.0e-3   1.0e-2  -2.0e-3   3.0e-3   0.0   0.0
  -1.5e-3  -4.0e-3   1.0e-3  -8.0e-3   0.0   1.0e-3
]

"""A parameterized custom field function returning physical fields."""
function uniform_magnetic_field(x, y, s, t, parameters)
  carrier = zero(x)
  return EMField(carrier, carrier, carrier,
                 carrier + parameters.Bx,
                 carrier + parameters.By,
                 carrier + parameters.Bz)
end

"""A callable field function; its fifth argument is supplied even without parameters."""
struct UniformBy{T}
  strength::T
end

function (em_field::UniformBy)(x, y, s, t, parameters)
  carrier = zero(x)
  return EMField(carrier, carrier, carrier,
                 carrier, carrier + em_field.strength, carrier)
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

  println("Runge-Kutta additive field-function demo")
  particle_count = size(INITIAL_COORDS, 1)
  println("Tracking $particle_count electrons at 1 GeV through $(length) m elements.\n")

  # 1. With no field contributions, RK reproduces exact drift tracking.
  exact_drift = track_element(Drift(L=length), species, p_over_q_ref)
  rk_drift = track_element(
    Drift(L=length, tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  verify_case("RK vs exact drift", rk_drift, exact_drift)

  # 2. A custom physical field agrees with the standard dipole parameters.
  by = 4.0e-3
  parameters = (Bx=0.0, By=by, Bz=0.0)
  functional_result = track_element(
    Drift(L=length, em_field=uniform_magnetic_field,
          em_field_params=parameters,
          tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  dipole_result = track_element(
    Drift(L=length, Kn0=by / p_over_q_ref,
          tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  verify_case("Field function vs standard dipole", functional_result, dipole_result)

  # 3. The whole-group API can declare rigidity-normalized custom fields.
  normalized_result = track_element(
    Drift(L=length,
          EMFieldParams=EMFieldParams(
            em_field=uniform_magnetic_field,
            em_field_params=(Bx=0.0, By=by / p_over_q_ref, Bz=0.0),
            em_field_normalized=true),
          tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  verify_case("Normalized vs physical function", normalized_result, functional_result)

  # 4. A custom dipole contribution adds to the standard quadrupole field.
  kn1 = 0.2
  combined_element = track_element(
    Quadrupole(L=length, Kn1=kn1,
               em_field=uniform_magnetic_field,
               em_field_params=parameters,
               tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  combined_multipoles = track_element(
    Quadrupole(L=length, Kn1=kn1, Kn0=by / p_over_q_ref,
               tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  verify_case("Additive function vs multipoles", combined_element, combined_multipoles)

  # 5. A callable object also receives the fifth argument (nothing by default).
  custom_result = track_element(
    Drift(L=length, em_field=UniformBy(by),
          tracking_method=RungeKutta(n_steps=n_steps)),
    species, p_over_q_ref)
  verify_case("Callable vs parameterized function", custom_result, functional_result)

  println("\nAll field-function cases passed.")
  @printf(
    "Combined-field final coordinates (particle 1):\n  x=% .6e  px=% .6e  y=% .6e  py=% .6e  zeta=% .6e  delta=% .6e\n",
    combined_element[1, :]...,
  )
  return nothing
end

main()
