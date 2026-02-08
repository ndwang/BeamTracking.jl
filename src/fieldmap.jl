# =====================================================================
# Field Map Data Structures, Interpolation, and Field Evaluation
# =====================================================================

# ---- Grid Descriptor Types ----

"""
    RectGrid3D{T}

3D rectangular grid descriptor for field maps.
Stores openPMD-beamphysics grid attributes.
"""
struct RectGrid3D{T}
  gridOriginOffset::NTuple{3, T}   # Physical coordinates of grid origin
  gridSpacing::NTuple{3, T}        # (dx, dy, dz)
  gridLowerBound::NTuple{3, Int}   # Lower bound index per dimension
  gridSize::NTuple{3, Int}         # (nx, ny, nz)
end

"""
    CylGrid2D{T}

2D cylindrical (r, z) grid descriptor for azimuthally symmetric field maps.
Stores openPMD-beamphysics grid attributes.
"""
struct CylGrid2D{T}
  gridOriginOffset::NTuple{2, T}   # (r_origin, z_origin)
  gridSpacing::NTuple{2, T}        # (dr, dz)
  gridLowerBound::NTuple{2, Int}   # Lower bound indices
  gridSize::NTuple{2, Int}         # (nr, nz)
end

"""
    grid_min(g, i)

Compute the minimum physical coordinate along dimension `i` for grid `g`.
"""
@inline grid_min(g, i) = g.gridOriginOffset[i] + g.gridLowerBound[i] * g.gridSpacing[i]

# ---- FieldMap Struct ----

"""
    FieldMap{Grid, E_re, E_im, B_re, B_im}

Field map for RK4 particle tracking through gridded electromagnetic fields.

Type parameters are each `Nothing` or a concrete array type, resolved at compile time
for zero-cost dispatch. Static fields have `E_im=nothing`, `B_im=nothing`.

Array layout: component-first — `(3, nx, ny, nz)` for 3D rectangular,
`(3, nr, nz)` for 2D cylindrical. All 3 vector components are contiguous at each grid point.
"""
struct FieldMap{Grid, E_re, E_im, B_re, B_im}
  grid::Grid
  E_re::E_re           # Nothing or AbstractArray — shape (3, spatial_dims...)
  E_im::E_im           # Nothing or AbstractArray (RF imaginary part)
  B_re::B_re           # Nothing or AbstractArray
  B_im::B_im           # Nothing or AbstractArray
  frequency::Float64   # ω (rad/s), 0 for static
  phase::Float64       # RF phase (rad)
  eleAnchorPt::Symbol  # :beginning, :center, :end (openPMD attribute)
  harmonic::Int        # 0 for static, >0 for RF (openPMD attribute)
end

Adapt.@adapt_structure FieldMap

# ---- Constructors ----

"""
    FieldMap(grid::RectGrid3D; E=nothing, B=nothing, frequency=0.0, phase=0.0,
             eleAnchorPt=:beginning, harmonic=0)

Construct a static rectangular field map. At least one of `E` or `B` must be provided.
Arrays must have shape `(3, nx, ny, nz)`.
"""
function FieldMap(grid::RectGrid3D; E=nothing, B=nothing, frequency=0.0, phase=0.0,
                  eleAnchorPt::Symbol=:beginning, harmonic::Int=0)
  nx, ny, nz = grid.gridSize
  E !== nothing && @assert size(E) == (3, nx, ny, nz) "E shape must be (3,$nx,$ny,$nz), got $(size(E))"
  B !== nothing && @assert size(B) == (3, nx, ny, nz) "B shape must be (3,$nx,$ny,$nz), got $(size(B))"
  E === nothing && B === nothing && error("At least one of E or B must be provided")
  return FieldMap(grid, E, nothing, B, nothing, Float64(frequency), Float64(phase), eleAnchorPt, harmonic)
end

"""
    FieldMap(grid::CylGrid2D; E=nothing, B=nothing, frequency=0.0, phase=0.0,
             eleAnchorPt=:beginning, harmonic=0)

Construct a static cylindrical field map. At least one of `E` or `B` must be provided.
Arrays must have shape `(3, nr, nz)`.
"""
function FieldMap(grid::CylGrid2D; E=nothing, B=nothing, frequency=0.0, phase=0.0,
                  eleAnchorPt::Symbol=:beginning, harmonic::Int=0)
  nr, nz = grid.gridSize
  E !== nothing && @assert size(E) == (3, nr, nz) "E shape must be (3,$nr,$nz), got $(size(E))"
  B !== nothing && @assert size(B) == (3, nr, nz) "B shape must be (3,$nr,$nz), got $(size(B))"
  E === nothing && B === nothing && error("At least one of E or B must be provided")
  return FieldMap(grid, E, nothing, B, nothing, Float64(frequency), Float64(phase), eleAnchorPt, harmonic)
end

# ---- Trilinear Interpolation (3D Rectangular) ----

"""
    interp_component(c, arr, grid::RectGrid3D, x, y, z)

Interpolate component `c` (1, 2, or 3) from a `(3, nx, ny, nz)` array using
trilinear interpolation. Branchless, returns zero for out-of-bounds points.
"""
@inline function interp_component(c::Int, arr, grid::RectGrid3D, x, y, z)
  x_min = grid_min(grid, 1);  dx = grid.gridSpacing[1];  nx = grid.gridSize[1]
  y_min = grid_min(grid, 2);  dy = grid.gridSpacing[2];  ny = grid.gridSize[2]
  z_min = grid_min(grid, 3);  dz = grid.gridSpacing[3];  nz = grid.gridSize[3]

  # Fractional 0-based indices
  fx = (x - x_min) / dx
  fy = (y - y_min) / dy
  fz = (z - z_min) / dz

  # Bounds check (branchless)
  in_bounds = (fx >= 0) & (fx <= nx - 1) &
              (fy >= 0) & (fy <= ny - 1) &
              (fz >= 0) & (fz <= nz - 1)

  # Clamp fractional index to [0, n-1], then clamp cell index to ensure ix+1 <= n
  fx_s = clamp(fx, zero(fx), oftype(fx, nx - 1))
  fy_s = clamp(fy, zero(fy), oftype(fy, ny - 1))
  fz_s = clamp(fz, zero(fz), oftype(fz, nz - 1))

  ix = min(unsafe_trunc(Int, fx_s) + 1, nx - 1)
  iy = min(unsafe_trunc(Int, fy_s) + 1, ny - 1)
  iz = min(unsafe_trunc(Int, fz_s) + 1, nz - 1)
  tx = fx_s - (ix - 1);  ty = fy_s - (iy - 1);  tz = fz_s - (iz - 1)

  # Trilinear — 8 corners in (3, nx, ny, nz) layout
  val = arr[c,ix,  iy,  iz  ]*(1-tx)*(1-ty)*(1-tz) +
        arr[c,ix+1,iy,  iz  ]*tx*(1-ty)*(1-tz) +
        arr[c,ix,  iy+1,iz  ]*(1-tx)*ty*(1-tz) +
        arr[c,ix+1,iy+1,iz  ]*tx*ty*(1-tz) +
        arr[c,ix,  iy,  iz+1]*(1-tx)*(1-ty)*tz +
        arr[c,ix+1,iy,  iz+1]*tx*(1-ty)*tz +
        arr[c,ix,  iy+1,iz+1]*(1-tx)*ty*tz +
        arr[c,ix+1,iy+1,iz+1]*tx*ty*tz

  return vifelse(in_bounds, val, zero(val))
end

"""
    interp_field(arr, grid::RectGrid3D, x, y, z)

Interpolate all 3 components from a `(3, nx, ny, nz)` array using trilinear interpolation.
Returns a tuple `(F1, F2, F3)`.
"""
@inline function interp_field(arr, grid::RectGrid3D, x, y, z)
  return (interp_component(1, arr, grid, x, y, z),
          interp_component(2, arr, grid, x, y, z),
          interp_component(3, arr, grid, x, y, z))
end

# ---- Bilinear Interpolation (2D Cylindrical) ----

"""
    interp_component(c, arr, grid::CylGrid2D, r, z)

Interpolate component `c` (1, 2, or 3) from a `(3, nr, nz)` array using
bilinear interpolation. Branchless, returns zero for out-of-bounds points.
"""
@inline function interp_component(c::Int, arr, grid::CylGrid2D, r, z)
  r_min = grid_min(grid, 1);  dr = grid.gridSpacing[1];  nr = grid.gridSize[1]
  z_min = grid_min(grid, 2);  dz = grid.gridSpacing[2];  nz = grid.gridSize[2]

  # Fractional 0-based indices
  fr = (r - r_min) / dr
  fz = (z - z_min) / dz

  # Bounds check (branchless)
  in_bounds = (fr >= 0) & (fr <= nr - 1) &
              (fz >= 0) & (fz <= nz - 1)

  # Clamp fractional index to [0, n-1], then clamp cell index to ensure ir+1 <= n
  fr_s = clamp(fr, zero(fr), oftype(fr, nr - 1))
  fz_s = clamp(fz, zero(fz), oftype(fz, nz - 1))

  ir = min(unsafe_trunc(Int, fr_s) + 1, nr - 1)
  iz = min(unsafe_trunc(Int, fz_s) + 1, nz - 1)
  tr = fr_s - (ir - 1);  tz = fz_s - (iz - 1)

  # Bilinear — 4 corners in (3, nr, nz) layout
  val = arr[c,ir,  iz  ]*(1-tr)*(1-tz) +
        arr[c,ir+1,iz  ]*tr*(1-tz) +
        arr[c,ir,  iz+1]*(1-tr)*tz +
        arr[c,ir+1,iz+1]*tr*tz

  return vifelse(in_bounds, val, zero(val))
end

"""
    interp_field(arr, grid::CylGrid2D, r, z)

Interpolate all 3 components from a `(3, nr, nz)` array using bilinear interpolation.
Returns a tuple `(F1, F2, F3)`.
"""
@inline function interp_field(arr, grid::CylGrid2D, r, z)
  return (interp_component(1, arr, grid, r, z),
          interp_component(2, arr, grid, r, z),
          interp_component(3, arr, grid, r, z))
end

# ---- Nothing-aware interpolation helpers (compile-time dispatch) ----

@inline _interp_field(::Nothing, grid::RectGrid3D, x, y, z) = (zero(x), zero(x), zero(x))
@inline _interp_field(arr, grid::RectGrid3D, x, y, z) = interp_field(arr, grid, x, y, z)

@inline _interp_field(::Nothing, grid::CylGrid2D, r, z) = (zero(r), zero(r), zero(r))
@inline _interp_field(arr, grid::CylGrid2D, r, z) = interp_field(arr, grid, r, z)

# ---- Field Evaluation: Rectangular (static) ----

"""
    fieldmap_em_field(x, y, z, pz, s, fm::FieldMap{<:RectGrid3D, ...}, z_offset)

Evaluate static rectangular field map at position `(x, y, s + z_offset)`.
Returns `(Ex, Ey, Ez, Bx, By, Bz)`.
"""
@inline function fieldmap_em_field(x, y, z, pz, s,
    fm::FieldMap{<:RectGrid3D, E_re, Nothing, B_re, Nothing},
    z_offset) where {E_re, B_re}
  z_map = s + z_offset
  Ex, Ey, Ez = _interp_field(fm.E_re, fm.grid, x, y, z_map)
  Bx, By, Bz = _interp_field(fm.B_re, fm.grid, x, y, z_map)
  return (Ex, Ey, Ez, Bx, By, Bz)
end

# ---- Field Evaluation: Cylindrical (static) ----

"""
    fieldmap_em_field(x, y, z, pz, s, fm::FieldMap{<:CylGrid2D, ...}, z_offset)

Evaluate static cylindrical field map at position `(x, y, s + z_offset)`.
Interpolates in cylindrical `(r, z)` coordinates, then converts to Cartesian.
At `r=0`, `Fr=Fθ=0` by azimuthal symmetry, so `cos_t`/`sin_t` are irrelevant.
Returns `(Ex, Ey, Ez, Bx, By, Bz)`.
"""
@inline function fieldmap_em_field(x, y, z, pz, s,
    fm::FieldMap{<:CylGrid2D, E_re, Nothing, B_re, Nothing},
    z_offset) where {E_re, B_re}
  z_map = s + z_offset
  r = sqrt(x^2 + y^2)

  # Exact trig from coordinates. At r=0, Fr=Fθ=0 by azimuthal symmetry,
  # so cos_t/sin_t are irrelevant (multiplied by zero).
  r_safe = vifelse(r > zero(r), r, one(r))
  cos_t = x / r_safe
  sin_t = y / r_safe

  # Interpolate in cylindrical basis: (Fr, Fθ, Fz)
  Er, Et, Ez = _interp_field(fm.E_re, fm.grid, r, z_map)
  Br, Bt, Bz = _interp_field(fm.B_re, fm.grid, r, z_map)

  # Exact rotation to Cartesian basis
  Ex = Er * cos_t - Et * sin_t
  Ey = Er * sin_t + Et * cos_t
  Bx = Br * cos_t - Bt * sin_t
  By = Br * sin_t + Bt * cos_t
  return (Ex, Ey, Ez, Bx, By, Bz)
end

# ---- Field Evaluation: RF (stub) ----

# RF field map tracking (harmonic > 0) — not yet implemented.
# Pending coordination on RF time convention.
@inline function fieldmap_em_field(x, y, z, pz, s,
    fm::FieldMap{G, E_re, E_im, B_re, B_im}, z_offset) where {G, E_re, E_im, B_re, B_im}
  error("RF field map tracking not yet implemented")
end
