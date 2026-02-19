# Field Map Support

## Overview

The field map system allows RK4 tracking through gridded electromagnetic field data — fields imported from electromagnetic solvers or measurements rather than described by analytic multipole expansions. It supports 3D rectangular and 2D cylindrically symmetric grids, with conventions following the openPMD-beamphysics standard.

All code is in `src/fieldmap.jl` (282 lines), with tests in `test/fieldmap_test.jl` (658 lines).

---

## Grid Descriptors

Two grid types describe the spatial layout of the field data:

### `RectGrid3D{T}` — 3D Cartesian grid

```julia
struct RectGrid3D{T}
  gridOriginOffset::NTuple{3, T}   # (x₀, y₀, z₀) physical origin
  gridSpacing::NTuple{3, T}        # (dx, dy, dz)
  gridLowerBound::NTuple{3, Int}   # lower bound index per dimension
  gridSize::NTuple{3, Int}         # (nx, ny, nz) number of grid points
end
```

The minimum physical coordinate along dimension `i` is:

```
grid_min(g, i) = gridOriginOffset[i] + gridLowerBound[i] * gridSpacing[i]
```

### `CylGrid2D{T}` — 2D cylindrical grid (azimuthally symmetric)

```julia
struct CylGrid2D{T}
  gridOriginOffset::NTuple{2, T}   # (r₀, z₀)
  gridSpacing::NTuple{2, T}        # (dr, dz)
  gridLowerBound::NTuple{2, Int}
  gridSize::NTuple{2, Int}         # (nr, nz)
end
```

Same `grid_min` formula applies.

---

## The `FieldMap` Struct

```julia
struct FieldMap{Grid, E_re, E_im, B_re, B_im}
  grid::Grid
  E_re::E_re       # Nothing or array
  E_im::E_im       # Nothing or array (RF imaginary part)
  B_re::B_re       # Nothing or array
  B_im::B_im       # Nothing or array (RF imaginary part)
  frequency::Float64
  phase::Float64
  eleAnchorPt::Symbol   # :beginning, :center, or :end
  harmonic::Int         # 0 for static, >0 for RF
end
```

Key design choices:

- **Component-first array layout**: `(3, nx, ny, nz)` for 3D, `(3, nr, nz)` for 2D. All 3 vector components are contiguous at each grid point.
- **Absent fields are `Nothing`**: Each of `E_re`, `E_im`, `B_re`, `B_im` is a type parameter. When a component is `nothing`, the compiler eliminates its interpolation entirely — zero runtime cost for absent fields.
- **GPU-compatible**: Uses `Adapt.@adapt_structure` so the struct can be transferred to GPU memory.
- **Constructor validates shapes** and requires at least one of `E` or `B`.

### Construction examples

```julia
# 3D rectangular, B-only (static magnet)
grid = RectGrid3D((0.0, 0.0, 0.0), (0.001, 0.001, 0.01), (0, 0, 0), (50, 50, 100))
B = zeros(3, 50, 50, 100)  # fill with field data...
fm = FieldMap(grid; B=B)

# 2D cylindrical, E-only (RF cavity)
grid = CylGrid2D((0.0, 0.0), (0.001, 0.01), (0, 0), (30, 200))
E = zeros(3, 30, 200)  # fill with field data...
fm = FieldMap(grid; E=E)

# Both E and B, with anchor point
fm = FieldMap(grid; E=E, B=B, eleAnchorPt=:center)
```

---

## Interpolation

### Trilinear — 3D rectangular grids

`interp_component(c, arr, grid::RectGrid3D, x, y, z)` interpolates component `c` (1, 2, or 3):

1. **Fractional indices**: Convert physical `(x, y, z)` to 0-based fractional grid indices:
   ```
   fx = (x - x_min) / dx
   ```

2. **Branchless bounds check**: Uses bitwise `&` (not short-circuit `&&`) so it works under SIMD/GPU:
   ```julia
   in_bounds = (fx >= 0) & (fx <= nx-1) & (fy >= 0) & (fy <= ny-1) & (fz >= 0) & (fz <= nz-1)
   ```

3. **Index clamping**: Clamp fractional index to `[0, n-1]`, then compute the 1-based cell index:
   ```julia
   fx_s = clamp(fx, 0, nx-1)
   ix = min(unsafe_trunc(Int, fx_s) + 1, nx - 1)   # ensures ix+1 <= nx
   tx = fx_s - (ix - 1)                              # fractional part within cell
   ```
   The `min(..., nx-1)` handles the boundary grid point (where `fx == nx-1` exactly) by pulling `ix` back so `ix+1` stays in bounds.

4. **8-corner trilinear weighting**:
   ```julia
   val = arr[c,ix,  iy,  iz  ]*(1-tx)*(1-ty)*(1-tz) +
         arr[c,ix+1,iy,  iz  ]*tx*(1-ty)*(1-tz) +
         ...  # all 8 corners
   ```

5. **Zero for out-of-bounds** (branchless):
   ```julia
   return vifelse(in_bounds, val, zero(val))
   ```

`interp_field(arr, grid, x, y, z)` calls `interp_component` for all 3 components, returning a tuple `(F1, F2, F3)`.

### Bilinear — 2D cylindrical grids

Same algorithm with 4-corner bilinear weights in `(r, z)`. The dispatch is on `grid::CylGrid2D`.

### Nothing-aware helpers

When a field component is `nothing`, dispatch skips interpolation entirely:

```julia
_interp_field(::Nothing, grid::RectGrid3D, x, y, z) = (zero(x), zero(x), zero(x))
_interp_field(arr,       grid::RectGrid3D, x, y, z) = interp_field(arr, grid, x, y, z)
```

This is resolved at compile time since `Nothing` vs array is encoded in the type parameter.

---

## Field Evaluation

### Static rectangular

```julia
fieldmap_em_field(x, y, z, pz, s, fm::FieldMap{<:RectGrid3D, ...}, z_offset)
```

- Evaluates at map coordinate `z_map = s + z_offset`, where `z_offset` accounts for `eleAnchorPt` (e.g. `z_offset = -L/2` for `:center`)
- Interpolates `E_re` and `B_re` via `_interp_field` (absent fields return zeros at no cost)
- Returns `(Ex, Ey, Ez, Bx, By, Bz)` in physical units (V/m, Tesla)

### Static cylindrical

```julia
fieldmap_em_field(x, y, z, pz, s, fm::FieldMap{<:CylGrid2D, ...}, z_offset)
```

The cylindrical case has an extra step — coordinate conversion:

1. Compute `r = sqrt(x² + y²)`
2. Handle `r = 0` safely: `r_safe = vifelse(r > 0, r, 1)` — at `r = 0`, `Fr` and `Fθ` are zero by azimuthal symmetry, so `cos_t`/`sin_t` multiply zero and their values don't matter
3. Interpolate in cylindrical basis: `(Fr, Fθ, Fz)` from `_interp_field(arr, grid, r, z_map)`
4. Rotate to Cartesian:
   ```
   Fx = Fr·cos(θ) − Fθ·sin(θ)
   Fy = Fr·sin(θ) + Fθ·cos(θ)
   ```
   where `cos(θ) = x/r_safe`, `sin(θ) = y/r_safe` (exact trig from coordinates, no `atan` call)

### RF fields

A stub method exists for field maps with nonzero imaginary parts (`E_im` or `B_im`), but it currently raises an error. This is pending coordination on the RF time convention.

---

## `FieldMapSource` — Bridging Field Maps to the RK4 Stepper

### Purpose

The RK4 kernel (`rk4_step!`, `rk4_kernel!`) is written generically over a "field source" — it calls `eval_em_field(field, x, y, z, pz, s)` without knowing where the fields come from. `FieldMapSource` is the adapter that makes a `FieldMap` satisfy this interface.

Its counterpart is `MultipoleSource`, which evaluates fields from analytic multipole coefficients. Both implement the same two-method interface, so the same RK4 code handles both.

### Definition

```julia
struct FieldMapSource{FM, T}
  fieldmap::FM        # FieldMap{Grid, E_re, E_im, B_re, B_im}
  z_offset::T         # longitudinal offset (from eleAnchorPt)
  g_bend::T           # curvature 1/ρ (0.0 for straight elements)
end
```

**`fieldmap`**: The `FieldMap` instance containing the grid descriptor and field arrays. Its full type (grid geometry, which field components are present) is encoded in the `FM` type parameter, so all dispatch decisions are resolved at compile time.

**`z_offset`**: A constant shift applied to the longitudinal coordinate before field lookup. During tracking, the RK4 kernel steps `s` from `0` to `L`. The field map may be defined in a different frame depending on `eleAnchorPt`:

| `eleAnchorPt` | `z_offset` | `z_map = s + z_offset` range |
|---|---|---|
| `:beginning` | `0` | `[0, L]` — map z=0 at element entrance |
| `:center` | `-L/2` | `[-L/2, L/2]` — map z=0 at element center |
| `:end` | `-L` | `[-L, 0]` — map z=0 at element exit |

This is computed once when constructing `FieldMapSource`, not during tracking. Currently, the caller computes it manually; a future Beamlines.jl integration would read `fm.eleAnchorPt` and compute it automatically.

**`g_bend`**: Element curvature for the equations of motion. For field-map elements this is typically `0.0` (straight geometry), but the field is included for interface consistency with `MultipoleSource`.

### Interface

`FieldMapSource` implements the two methods the RK4 stepper requires:

```julia
# Field evaluation — delegates to fieldmap_em_field with the stored z_offset
@inline eval_em_field(f::FieldMapSource, x, y, z, pz, s) =
    fieldmap_em_field(x, y, z, pz, s, f.fieldmap, f.z_offset)

# Curvature — returns the stored g_bend
@inline get_g_bend(f::FieldMapSource) = f.g_bend
```

Because `FieldMapSource` is a type parameter of `rk4_kernel!`, the compiler specializes the entire RK4 loop for the specific field map configuration (grid type, which E/B components exist). There is no dynamic dispatch and no overhead beyond the interpolation arithmetic itself.

### Usage

```julia
# Build grid and field map
grid = RectGrid3D((-0.02, -0.02, 0.0), (0.01, 0.01, 0.1), (0,0,0), (5, 5, 20))
B = zeros(3, 5, 5, 20)
B[3, :, :, :] .= 0.01  # uniform 0.01 T solenoid
fm = FieldMap(grid; B=B, eleAnchorPt=:beginning)

# Wrap in FieldMapSource with z_offset computed from eleAnchorPt
z_offset = 0.0  # :beginning → no shift
field = FieldMapSource(fm, z_offset, 0.0)

# Pass to the RK4 kernel — same call signature as MultipoleSource
RungeKuttaTracking.rk4_kernel!(i, coords, beta_0, tilde_m,
    charge, p0c, mc2, (0.0, L), ds_step, field)
```

### Comparison with `MultipoleSource`

| | `MultipoleSource` | `FieldMapSource` |
|---|---|---|
| Field origin | Analytic multipole coefficients | Gridded arrays + interpolation |
| Parameters | `mm`, `kn`, `ks`, `p_over_q_ref` | `FieldMap`, `z_offset` |
| Field units | Normalized → multiplied by `p_over_q_ref` | Already in physical units (T, V/m) |
| Typical use | Standard magnets (quad, sextupole, bend) | Imported solver/measurement data |
| `g_bend` | From element's `BendParams` | Typically `0.0` |

Both are concrete structs with fully parameterized types, so the compiler generates separate specialized code paths for each — no runtime branching between them.

---

## Data Flow Summary

```
Physical (x, y, s) coordinates
        │
        ▼
┌─────────────────────┐
│  fieldmap_em_field   │  ← dispatches on grid type (RectGrid3D or CylGrid2D)
│                      │     and field presence (Nothing vs array)
│  1. z_map = s + off  │
│  2. (cyl only)       │
│     r = √(x²+y²)    │
│  3. interp E, B      │  ← trilinear or bilinear, branchless
│  4. (cyl only)       │
│     rotate to Cart.  │
└──────────┬──────────┘
           │
    (Ex,Ey,Ez,Bx,By,Bz) in physical units
           │
           ▼
    kick_vector() → du/ds
           │
           ▼
    rk4_step!() → updated particle state
```

---

## Testing Highlights

The field map tests (`test/fieldmap_test.jl`) validate:

- **Interpolation accuracy**: linear functions are reproduced exactly, on-grid points recover stored values, out-of-bounds returns zero
- **Boundary handling**: single-cell grids, edge grid points
- **Type stability**: `@inferred` on all interpolation and field evaluation paths
- **Field map vs multipole agreement**: uniform `Bz` solenoid and quadrupole gradient tracked via field map match the same fields tracked via `MultipoleSource` (within interpolation error, `rtol ≈ 1e-3` to `1e-6`)
- **Physics validation**: energy gain in uniform `Ez`, transverse momentum conservation in uniform `Bz`, Larmor rotation in cylindrical solenoid
- **Cylindrical edge cases**: `r = 0` produces no NaN/Inf; radial field `Br = r` correctly maps to `Bx = x, By = y`
- **`eleAnchorPt`**: offset semantics are tested for `:beginning`, `:center`, `:end`
