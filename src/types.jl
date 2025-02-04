#=

Bunch, Coord, and Particle type definitions. StructArrays is 
used to handle an internal SoA layout of memory which also 
allows us to mutate, so both the phase space Coord struct 
(defined here) and Quaternion struct (defined in 
ReferenceFrameRotations) can be static. 

Particle struct simply goes from StructArrays SoA to AoS.

=#


"""
    Coord{T} <: FieldVector{6, T}

A static 6D phase space coordinate vector representing a particle's position and momentum.

# Fields
- `x::T`: Horizontal position
- `px::T`: Horizontal momentum normalized by reference beta*gamma
- `y::T`: Vertical position  
- `py::T`: Vertical momentum normalized by reference beta*gamma
- `z::T`: Longitudinal position
- `pz::T`: Longitudinal momentum normalized by reference beta*gamma

All fields default to 0.0 if not specified.
"""
Base.@kwdef struct Coord{T} <: FieldVector{6, T} 
  x::T  = 0.0  # Horizontal position
  px::T = 0.0  # Horizontal momentum 
  y::T  = 0.0  # Vertical position
  py::T = 0.0  # Vertical momentum
  z::T  = 0.0  # Longitudinal position
  pz::T = 0.0  # Longitudinal momentum
end

# Static quaternion type defined by ReferenceFrameRotatio

"""
    Bunch{T<:StructVector{<:Coord}, U<:Union{Nothing, StructVector{<:Quaternion}}}

A collection of particles with their phase space coordinates and optional spin quaternions.
Uses StructArrays for efficient memory layout (Structure of Arrays).

# Fields
- `species::Species`: Particle species (e.g., electron, proton)
- `beta_gamma_ref::Float64`: Reference β*γ value used to normalize momenta
- `v::T`: StructVector of phase space coordinates for all particles
- `q::U`: Optional StructVector of spin quaternions (nothing if spin tracking disabled)
"""
struct Bunch{T<:StructVector{<:Coord}, U<:Union{Nothing, StructVector{<:Quaternion}}}
  species::Species
  beta_gamma_ref::Float64
  v::T
  q::U
end

# Initialize a Bunch with either a single particle (scalars)
"""
    Bunch(; species::Species=Species("electron"), beta_gamma_ref=1.0, 
           spin::Union{Bool,Nothing}=nothing, gtpsa_map::Union{Bool,Nothing}=nothing,
           x::Union{Number,AbstractVector}=0.0, px::Union{Number,AbstractVector}=0.0, 
           y::Union{Number,AbstractVector}=0.0, py::Union{Number,AbstractVector}=0.0, 
           z::Union{Number,AbstractVector}=0.0, pz::Union{Number,AbstractVector}=0.0 )


Initializes a `Bunch`. Any of the specified phase space coordinates may be scalar `Number`s or 
`Vector`s to store as a structure-of-arrays. Internally, all phase space coordinates are stored 
as `Vector`s. If all phase space coordinates are scalar `Number`s, then a `Bunch` is created with 
a single particle. If any of the coordinates are specified as `Vector`s, then all other scalar-
specified quantities are `fill`-ed as `Vector`s. For example, `Bunch(x=1.0, y=[2,3])` creates a 
bunch with two particles having the phase space coordinates `[1.0, 0.0, 2.0, 0.0, 0.0, 0.0]` 
and `[1.0, 0.0, 3.0, 0.0, 0.0, 0.0]`.

## Other keyword arguments
• `species`         -- Particle species, default is electron
• `beta_gamma_ref`  -- Reference Lorentz beta*gammma to normalize the momenta to
• `spin`            -- If true, spin tracking is turned on and a quaternion for each particle is tracked
• `gtpsa_map`       -- If true, GTPSA map tracking is used for each particle using the Descriptor defined in GTPSA.desc_current
"""
function Bunch(; species::Species=Species("electron"), beta_gamma_ref=1.0, 
                spin::Union{Bool,Nothing}=nothing, gtpsa_map::Union{Bool,Nothing}=nothing,
                x::Union{Number,AbstractVector}=0.0, px::Union{Number,AbstractVector}=0.0, 
                y::Union{Number,AbstractVector}=0.0, py::Union{Number,AbstractVector}=0.0, 
                z::Union{Number,AbstractVector}=0.0, pz::Union{Number,AbstractVector}=0.0 )
                
  idx_vector = findfirst(t->t isa AbstractVector, (x, px, y, py, z, pz))
  if isnothing(idx_vector)
    N_particle = 1
  else
    N_particle = length(getindex((x, px, y, py, z, pz), idx_vector))
  end

  T1 = promote_type(eltype(x), eltype(px), eltype(y),eltype(py), eltype(z), eltype(pz)) 
  if !isnothing(gtpsa_map)
    if gtpsa_map == true
      GTPSA.numvars(GTPSA.desc_current) == 6 || error("Invalid GTPSA Descriptor! Number of variables must be equal to 6.")
      T = promote_type(TPS64, T1)
    else
      error("For no GTPSA map tracking, please omit the gtpsa_map kwarg or set gtpsa_map=nothing. This is to ensure type stability.")
    end
  else
    T = T1
  end
  
  @inline function make_vec_T(T, vec_or_num, N_particle)
    if vec_or_num isa AbstractVector
      if eltype(vec_or_num) == T
        return vec_or_num
      else
        return T.(vec_or_num)
      end
    else
      if isimmutable(T)
        return fill(T(vec_or_num),  N_particle)
      else
        vec = Vector{T}(undef, N_particle)
        for i in eachindex(vec)
          vec[i] = T(vec_or_num)
        end
        return vec
      end
    end
  end

  x1  = make_vec_T(T, x,  N_particle)
  px1 = make_vec_T(T, px, N_particle)
  y1  = make_vec_T(T, y,  N_particle)
  py1 = make_vec_T(T, py, N_particle)
  z1  = make_vec_T(T, z,  N_particle)
  pz1 = make_vec_T(T, pz, N_particle)

  coords = (x1, px1, y1, py1, z1, pz1)
  if !isnothing(gtpsa_map) # Set slopes if GTPSA map tracking
    for var_idx in eachindex(coords)
      for i in eachindex(coords[var_idx])
        coords[var_idx][i][var_idx] = 1.0
      end
    end
  end

  v = StructArray{Coord{T}}((x1, px1, y1, py1, z1, pz1))

  if !isnothing(spin)
    if spin == true
      q0 = make_vec_T(T, 1, N_particle)
      q1 = make_vec_T(T, 0, N_particle)
      q2 = make_vec_T(T, 0, N_particle)
      q3 = make_vec_T(T, 0, N_particle)
      q = StructArray{Quaternion{T}}((q0, q1, q2, q3))
    else
      error("For no spin tracking, please omit the spin kwarg or set spin=nothing. This is to ensure type stability.")
    end
  else
    q = nothing
  end


  return Bunch(species, Float64(beta_gamma_ref), v, q)
end

"""
    Particle{T,U<:Union{Nothing,Quaternion{T}}}

A single particle representation extracted from a Bunch.
Converts from Structure of Arrays (SoA) to Array of Structures (AoS) format.

# Fields
- `species::Species`: Particle species (e.g., electron, proton)
- `beta_gamma_ref::Float64`: Reference β*γ value used to normalize momenta
- `v::Coord{T}`: Phase space coordinates of the particle
- `q::U`: Optional spin quaternion (nothing if spin tracking disabled)
"""
struct Particle{T,U<:Union{Nothing,Quaternion{T}}}
  species::Species
  beta_gamma_ref::Float64
  v::Coord{T}
  q::U
end

"""
    Particle(bunch::Bunch, idx::Integer=1)

Construct a single Particle from a Bunch at the specified index (defaults to first particle).
Extracts the particle's coordinates and spin (if present) from the Bunch's StructArrays.

# Arguments
- `bunch::Bunch`: The bunch to extract the particle from
- `idx::Integer=1`: Index of the particle to extract (default: 1)

# Returns
- `Particle`: A new Particle instance with the extracted data
"""
function Particle(bunch::Bunch, idx::Integer=1)
  v = bunch.v[idx] # StructArrays handles this!
  q = isnothing(bunch.q) ? nothing : bunch.q[idx]
  return Particle(bunch.species, bunch.beta_gamma_ref, v, q)
end
