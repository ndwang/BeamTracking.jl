module Linear
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

"""
    struct Linear.Drift{T}
  
## Fields
• `L` -- Drift length / m \\
"""
Base.@kwdef struct Drift{T}
  L::T  # drift length / m
end

"""
    struct Linear.Quadrupole{T}

## Fields
• `L`  -- Quadrupole length / m \\
• `Bn1` -- Quadrupole gradient / (T·m^-1)
"""
Base.@kwdef struct Quadrupole{T}
  L::T   # quadrupole length / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end


Base.@kwdef struct SBend{T}
  L::T   # arc length / m
  B0::T  # field strength / T
  g::T   # coordinate system curvature through element / m^-1
  e1::T  # edge 1 angle / rad
  e2::T  # edge 2 angle / rad
end

Base.@kwdef struct Combined{T}
  L::T 
  B0::T 
  Bn1::T
  g::T
  e1::T
  e2::T
end

Base.@kwdef struct Solenoid{T}
  L::T
  Bs::T
end


function track!(bunch::Bunch, ele::Linear.Drift; work=get_work(bunch, Val{0}()))
  L = ele.L
  v = bunch.v
  
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)

  @FastGTPSA! begin
  @. v.x  = v.x + L*v.px
  @. v.y  = v.y + L*v.py
  @. v.z  = v.z + L/gamma_ref^2*v.pz
  end

  # Spin unchanged
  
  return bunch
end


function track!(bunch::Bunch, ele::Linear.Quadrupole; work=get_work(bunch, Val{1}()))
  v = bunch.v
  L = ele.L

  Kn1 = ele.Bn1 / brho(massof(bunch.species), bunch.beta_gamma_ref, chargeof(bunch.species))
  gamma_ref = sr_gamma(bunch.beta_gamma_ref)

  if Kn1 >= 0
    k = sqrt(Kn1)
    kL = k*L
    sx  = sin(kL)
    cx  = cos(kL)
    sxc = sincu(kL)
    sy  = sinh(kL)
    cy  = cosh(kL)
    syc = sinhcu(kL)
    sgn = 1
  else
    k = sqrt(-Kn1)
    kL = k*L
    sx  = sinh(kL)
    cx  = cosh(kL)
    sxc = sinhcu(kL)
    sy  = sin(kL)
    cy  = cos(kL)
    syc = sincu(kL)
    sgn = -1
  end

  
  # Note: 0+work[1] is workaround for GTPSA @FastGTPSA! bug
  @FastGTPSA! begin
  @. work[1]  = cx*v.x + sxc*L*v.px
  @. v.px     = -sgn*k*sx*v.x + cx*v.px
  @. v.x      = 0+work[1] 

  @. work[1]  = cy*v.y + syc*L*v.py
  @. v.py     = sgn*k*sy*v.y + cy*v.py
  @. v.y      = 0+work[1] 

  @. v.z      = v.z + L/gamma_ref^2*v.pz
  end 

  return bunch
end 

"""
Routine to linearly tracking through a Sector Magnet
"""
function track!(beamf::Bunch, ele::Linear.SBend, beami::Bunch)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.vec
  zf = beamf.vec
  L = ele.L
  q = chargeof(beami.species)

  #curvature of B field
  k = ele.B0 / brho(massof(beami.species),beami.beta_gamma_ref,q)
  kl = k * L


  @FastGTPSA! begin
  @. zf.x  = zi.x * cos(kl) +  zi.px * sin(kl) / k + zi.pz * (1 - cos(kl)) / k
  @. zf.px = -zi.x * k * sin(kl) + zi.px * cos(kl) + zi.pz * sin(kl)
  @. zf.y = zi.y + zi.py * L
  @. zf.py = zi.py
  @. zf.z = -zi.x * sin(kl) + zi.px * (cos(kl) - 1) / k + zi.z + zi.pz * (sin(kl) - kl) / k
  @. zf.pz = zi.pz
  end 
end 


"""
Routine to linearly tracking through a Combined Magnet
"""
function track!(beamf::Bunch, ele::Linear.Combined, beami::Bunch)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.vec
  zf = beamf.vec
  L = ele.L
  
  brho = brho(massof(beami.species),beami.beta_gamma_ref,q)
  ka = ele.B0 / brho #curvature
  k1 =  ele.Bn1 / brho #quad strength

  K = ka^2 + k1

  sqk1 = sqrt(abs(k1)) 
  sqK = sqrt(abs(K))

  Kl = sqK * L
  kl  = sqk1 * L

  @FastGTPSA! begin
  @. zf.x  = zi.x * cos(Kl) +  zi.px * sin(Kl) / sqK + zi.pz * (1 - cos(Kl)) * ka / K
  @. zf.px = - zi.x * sqK * sin(Kl) + zi.px * cos(Kl) + zi.pz * sin(Kl) * ka / sqK
  @. zf.y =  zi.y * cosh(kl) + zi.py * sinh(kl) / sqk
  @. zf.py = zi.y * sqK * sinh(kl) + zi.py * cosh(kl)
  @. zf.z = - zi.x * sin(Kl) * ka / sqK + zi.px * (cos(Kl)-1) * ka / K + zi.z + zi.pz * (sin(Kl) / sqK - l) *  ka^2 / K
  @. zf.pz = zi.pz
  end 
end 

end