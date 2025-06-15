
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]

# Thin elements
@inline thin_bsolenoid(tm::Linear, bunch, bm0) = error("Undefined for tracking method $tm")   
@inline thin_bdipole(tm::Linear, bunch, bm1)   = error("Undefined for tracking method $tm") 
@inline function thin_bquadrupole(tm::Linear, bunch, bm2)
  K1L = get_thin_strength(bm2, 0, bunch.Brho_ref)
  mx, my = LinearTracking.linear_thin_quad_matrices(K1L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, 0, nothing, nothing))
end
@inline thin_bmultipole(tm::Linear, bunch, bdict) = KernelCall() # Do nothing
@inline thin_bend_no_field(tm::Linear, bunch, bendparams)          = error("Undefined for tracking method $tm")                      
@inline thin_bend_bsolenoid(tm::Linear, bunch, bendparams, bm0)    = error("Undefined for tracking method $tm")                    
@inline thin_bend_bdipole(tm::Linear, bunch, bendparams, bm1)      = error("Undefined for tracking method $tm")                  
@inline thin_bend_bquadrupole(tm::Linear, bunch, bendparams, bm2)  = error("Undefined for tracking method $tm")                      
@inline thin_bend_bmultipole(tm::Linear, bunch, bendparams, bdict) = error("Undefined for tracking method $tm")                                           

# Thick elements
@inline function thick_bsolenoid(tm::Linear, bunch, bm0, L) 
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  Ks = get_thick_strength(bm0, L, bunch.Brho_ref)
  mxy = LinearTracking.linear_solenoid_matrix(Ks, L)
  return KernelCall(LinearTracking.linear_coast!, (mxy, L/gamma_0^2, nothing, nothing))
end

@inline thick_bdipole(tm::Linear, bunch, bm1, L)                       = error("Undefined for tracking method $tm")
@inline function thick_bquadrupole(tm::Linear, bunch, bm2, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  K1 = get_thick_strength(bm2, L, bunch.Brho_ref)
  mx, my = LinearTracking.linear_quad_matrices(K1, L)
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, L/gamma_0^2, nothing, nothing))
end

@inline function thick_bmultipole(tm::Linear, bunch, bdict, L)
  if length(keys(bdict)) > 1 
    error("Linear tracking currently does not support overlapping multipoles")
  else # higher order multipole - drift for these
    return drift(tm, bunch, L)
  end
end
@inline thick_bend_no_field(tm::Linear, bunch, bendparams, L)          = error("Undefined for tracking method $tm")                      
@inline thick_bend_bsolenoid(tm::Linear, bunch, bendparams, bm0, L)    = error("Undefined for tracking method $tm")                   
@inline function thick_bend_bdipole(tm::Linear, bunch, bendparams, bm1, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  K0 = get_thick_strength(bm1, L, bunch.Brho_ref)
  mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, gamma_0, bendparams.e1, bendparams.e2)
  if !(K0 ≈ bendparams.g)
    error("Linear tracking currently only supports BendParams.g ≈ BMultipoleParams.K0")
  end
  return KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t))
end
@inline thick_bend_bquadrupole(tm::Linear, bunch, bendparams, bm2, L)  = error("Undefined for tracking method $tm")                     
@inline thick_bend_bmultipole(tm::Linear, bunch, bendparams, bdict, L) = error("Undefined for tracking method $tm")                                          
@inline function drift(tm::Linear, bunch, L)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  return KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2))
end


#=
# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  b::BunchView,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::Linear;
  kwargs...
)
  # Unpack the line element (type unstable)
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L

  # Function barrier
  linear_universal!(i, b, bunch, L, bm, bp, ma; kwargs...)
end

@inline function get_thick_strength(bm, L, Brho_ref)
  s = bm.strength
  if !bm.normalized
    s /= Brho_ref
  end
  if bm.integrated
    if L == 0
      error("LineElement length is zero; cannot computed non-integrated strength")
    end
    s /= L
  end
  return s
end

@inline function get_thin_strength(bm, L, Brho_ref)
  s = bm.strength
  if !bm.normalized
    s /= Brho_ref
  end
  if !bm.integrated
    s *= L
  end
  return s
end

# Step 2: Push particles through -----------------------------------------
function linear_universal!(
  i, 
  b,
  bunch,
  L, 
  bmultipoleparams, 
  bendparams, 
  alignmentparams;
  kwargs...
) 
  kc = KernelChain(Val{1}())
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)
  if !isactive(bmultipoleparams) # Drift
    if isactive(bendparams)
      error("Linear tracking requires BendParams.g == BMultipoleParams.K0")
    end
    kc = push(kc, KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2)))
  elseif haskey(bmultipoleparams.bdict, 0) # Solenoid
    if any(t -> t >= 1, keys(bmultipoleparams.bdict))
      error("Linear tracking does not support combined solenoid + other multipole magnets")
    end
    if isactive(bendparams)
      error("Linear tracking does not currently support solenoid with bending")
    end
    if L == 0
      error("Thin solenoid not supported yet")
    end

    Ks = get_thick_strength(bmultipoleparams.bdict[0], L, bunch.Brho_ref)

    mxy = LinearTracking.linear_solenoid_matrix(Ks, L)
    kc = push(kc, KernelCall(LinearTracking.linear_coast!, (mxy, L/gamma_0^2, nothing, nothing)))
  elseif haskey(bmultipoleparams.bdict, 1) # Bend
    if !isactive(bendparams)
      error("Linear tracking requires BendParams.g ≈ BMultipoleParams.K0")
    end
    if L == 0
      error("Thin bend not supported yet")
    end
    if any(t -> t == 0 || t > 1, keys(bmultipoleparams.bdict))
      error("Combined function bend tracking not implemented yet")
    end

    K0 = get_thick_strength(bmultipoleparams.bdict[1], L, bunch.Brho_ref)

    if !(K0 ≈ bendparams.g)
      error("Linear tracking requires BendParams.g ≈ BMultipoleParams.K0")
    end
    mx, my, r56, d, t = LinearTracking.linear_bend_matrices(K0, L, gamma_0, bendparams.e1, bendparams.e2)
    kc = push(kc, KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56, d, t)))
  elseif haskey(bmultipoleparams.bdict, 2) # Quadrupole
    if isactive(bendparams)
      error("For Linear combined function magnet tracking, both the K0 multipole and BendParams must be set")
    end
    if L == 0
      K1L = get_thin_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
      mx, my = LinearTracking.linear_thin_quad_matrices(K1L)
    else
      K1 = get_thick_strength(bmultipoleparams.bdict[2], L, bunch.Brho_ref)
      mx, my = LinearTracking.linear_quad_matrices(K1, L)
    end
    kc = push(kc, KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, L/gamma_0^2, nothing, nothing)))
  else # Drift for higher-order multipoles
    kc = push(kc, KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2)))
  end

  runkernels!(i, b, kc; kwargs...)
  return nothing
end
=#