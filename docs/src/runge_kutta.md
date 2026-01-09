# Runge-Kutta Tracking

The `RungeKutta` tracking method provides a 4th-order Runge-Kutta (RK4) integration scheme for tracking particles through arbitrary electromagnetic fields. This method is particularly useful for elements with complex fields that cannot be handled analytically.

## Overview

The Runge-Kutta tracking implementation uses a classical 4th-order Runge-Kutta method to numerically integrate the relativistic equations of motion for charged particles in electromagnetic fields. The implementation is optimized for GPU/SIMD compatibility using branchless operations and stack-allocated StaticArrays.

## Configuration

The `RungeKutta` tracking method can be configured with the following parameters:

### Parameters

- **`ds_step`**: The step size for integration (in meters). If not specified, defaults to `0.2` meters when neither parameter is set.

- **`n_steps`** : The number of integration steps. If specified and positive, the step size is calculated as `h = L / n_steps` where `L` is the element length.

**Important**: Only one of `ds_step` or `n_steps` should be specified. If both are set to positive values, an error will be raised. If neither is specified, `ds_step=0.2` is used by default.

### Examples

```julia
# Use default step size (0.2 meters)
ele.tracking_method = RungeKutta()

# Specify step size explicitly
ele.tracking_method = RungeKutta(ds_step=0.1)

# Specify number of steps
ele.tracking_method = RungeKutta(n_steps=50)
```

## Physics

The Runge-Kutta method integrates the relativistic equations of motion expressed in terms of arc length `s` as the independent variable. The state vector is:

```math
\mathbf{u} = [x, p_x, y, p_y, z, p_z]
```

where:
- `x, y, z`: Position coordinates (meters)
- `p_x, p_y, p_z`: Normalized momentum components (relative to reference momentum `p₀`)

The derivatives `du/ds` are calculated from the relativistic electromagnetism:

### Position Derivatives

```math
\frac{dx}{ds} = v_x \frac{dt}{ds}, \quad \frac{dy}{ds} = v_y \frac{dt}{ds}
```

where `v_x, v_y` are the transverse velocity components and `dt/ds` accounts for the relationship between arc length and time.

### Momentum Derivatives

```math
\frac{dp_x}{ds} = \frac{F_x}{p₀c} \frac{dt}{ds} + g_{bend} \frac{p_z}{p₀}, \quad \frac{dp_y}{ds} = \frac{F_y}{p₀c} \frac{dt}{ds}
```

where `F_x, F_y` are the Lorentz force components and `g_bend` is the curvature parameter (nonzero only in bend elements).

### Longitudinal Coordinate

The `z` coordinate represents the phase space deviation and is integrated with corrections for energy changes and bending:

```math
\frac{dz}{ds} = \text{rel\_dir} \left(\frac{\beta}{\beta₀} - 1\right) + \text{rel\_dir} \frac{\sqrt{1-v_t²} - 1 - dh_{bend}}{\sqrt{1-v_t²}} + \frac{d\beta}{ds} \frac{z}{\beta}
```

### Energy Deviation

```math
\frac{dp_z}{ds} = \frac{dp}{ds} / p₀c
```

where `dp/ds` is calculated from the work done by the electromagnetic fields.

## Implementation Details

### RK4 Algorithm

The 4th-order Runge-Kutta method uses the standard four-stage scheme:

1. **k₁**: Evaluate derivatives at current state `u(s)`
2. **k₂**: Evaluate derivatives at `u(s) + (h/2) k₁`
3. **k₃**: Evaluate derivatives at `u(s) + (h/2) k₂`
4. **k₄**: Evaluate derivatives at `u(s) + h k₃`

The final update is:

```math
u(s+h) = u(s) + \frac{h}{6}(k₁ + 2k₂ + 2k₃ + k₄)
```

### Field Function Interface

The Runge-Kutta method requires a field function that returns the electromagnetic field components at a given position:

```julia
field_func(x, px, y, py, z, pz, s, field_params) -> (Ex, Ey, Ez, Bx, By, Bz)
```

When used with `Beamlines.jl` elements, the field function is automatically obtained from `Beamlines.field_calc(ele)`.

### Particle Loss Detection

The implementation includes automatic detection of unphysical particle states:

- **Unphysical momenta**: If the transverse velocity squared `v_t² ≥ 1`, the particle is marked as lost (`STATE_LOST_PZ`)
- The check is performed branchlessly for GPU compatibility
- Lost particles have their derivatives set to zero

### GPU/SIMD Compatibility

The implementation uses several techniques for GPU and SIMD compatibility:

- **Branchless operations**: Uses `vifelse` for conditional logic
- **Stack-allocated arrays**: All intermediate values use `SVector` from StaticArrays.jl
- **No dynamic memory allocation**: All arrays are stack-allocated
- **Kernel abstraction**: Compatible with `@makekernel` macro for GPU execution

## Usage

### Basic Usage with Beamlines.jl

```julia
using BeamTracking, Beamlines

# Create an element with Runge-Kutta tracking
ele = Quadrupole(L=1.0, k1=0.1)
ele.tracking_method = RungeKutta(ds_step=0.1)

# Create a bunch and track
bunch = Bunch(zeros(100, 6), R_ref=1e6, species=Species("electron"))
track!(bunch, ele)
```

### Custom Field Functions

For custom field configurations, you can use the low-level API:

```julia
using BeamTracking.RungeKuttaTracking

# Define a custom field function
function my_field(x, px, y, py, z, pz, s, params)
    Ex = params.E0 * sin(2π * s / params.λ)
    return (Ex, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# Setup parameters
species = Species("electron")
R_ref = 1e6
beta_0, gamsqr_0, tilde_m = drift_params(species, R_ref)
charge = chargeof(species) / E_CHARGE
p0c = R_to_pc(species, R_ref)
mc2 = massof(species)

# Track
s_span = (0.0, 1.0)
n_steps = 100
g_bend = 0.0
field_params = (E0=1e4, λ=0.1)

RungeKuttaTracking.rk4_kernel!(
    1, bunch.coords, beta_0, gamsqr_0, tilde_m,
    charge, p0c, mc2, s_span, n_steps, g_bend,
    my_field, field_params
)
```

## Supported Elements

The Runge-Kutta method works with all thick elements that have field definitions.
