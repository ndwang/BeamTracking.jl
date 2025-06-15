# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  b::BunchView,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  tm;
  kwargs...
)
  # Unpack the line element (type unstable)
  L = ele.L
  ap = ele.AlignmentParams
  bp = ele.BendParams
  bm = ele.BMultipoleParams
  pp = ele.PatchParams

  # Function barrier
  universal!(i, b, tm, bunch, L, ap, bp, bm, pp; kwargs...)
end

# Step 2: Push particles through -----------------------------------------
function universal!(
  i, 
  b,
  tm,
  bunch,
  L, 
  alignmentparams, 
  bendparams,
  bmultipoleparams,
  patchparams;
  kwargs...
) 
  kc = KernelChain(Val{1}())
  if isactive(alignmentparams)
    if isactive(patchparams)
      error("Tracking through a LineElement containing both PatchParams and AlignmentParams is undefined")
    end
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, true)))
  end

  if isactive(patchparams) 
    # Patch
    if isactive(bendparams) || isactive(bmultipoleparams)
      error("Tracking through a LineElement containing both PatchParams and BendParams/BMultipoleParams not currently defined")
    else
      kc = push(kc, @inline(pure_patch(tm, bunch, patchparams, L)))
    end
  elseif isactive(bendparams)
    # Bend
    if !isactive(bmultipoleparams) 
      # Bend no field
      kc = push(kc, @inline(bend_no_field(tm, bunch, bendparams, L)))
    else
      bdict = bmultipoleparams.bdict
      n_multipoles = length(keys(bdict))
      # BMultipole bend
      if n_multipoles > 1
        kc = push(kc, @inline(bend_bmultipole(tm, bunch, bendparams, bdict, L)))
      elseif haskey(bdict, 0)
        kc = push(kc, @inline(bend_bsolenoid(tm, bunch, bendparams, bdict[0], L)))
      elseif haskey(bdict, 1)
        kc = push(kc, @inline(bend_bdipole(tm, bunch, bendparams, bdict[1], L)))
      elseif haskey(bdict, 2)
        kc = push(kc, @inline(bend_bquadrupole(tm, bunch, bendparams, bdict[2], L)))
      else
        kc = push(kc, @inline(bend_bmultipole(tm, bunch, bendparams, bdict, L)))
      end
    end
  elseif isactive(bmultipoleparams)
    # BMultipole
    bdict = bmultipoleparams.bdict
    n_multipoles = length(keys(bdict))
    # BMultipole bend
    if n_multipoles > 1
      kc = push(kc, @inline(bmultipole(tm, bunch, bdict, L)))
    elseif haskey(bdict, 0)
      kc = push(kc, @inline(bsolenoid(tm, bunch, bdict[0], L)))
    elseif haskey(bdict, 1)
      kc = push(kc, @inline(bdipole(tm, bunch, bdict[1], L)))
    elseif haskey(bdict, 2)
      kc = push(kc, @inline(bquadrupole(tm, bunch, bdict[2], L)))
    else
      kc = push(kc, @inline(bmultipole(tm, bunch, bdict, L)))
    end
  elseif !(L ≈ 0)
    kc = push(kc, @inline(drift(tm, bunch, L)))
  end
    
  if isactive(alignmentparams)
    kc = push(kc, @inline(misalign(tm, bunch, alignmentparams, false)))
  end

  runkernels!(i, b, kc; kwargs...)
  return nothing
end

# Coordinate system transformations:
@inline misalign(tm, bunch, alignmentparams, in) = error("Undefined for tracking method $tm")
@inline pure_patch(tm, bunch, patchparams, L)    = error("Undefined for tracking method $tm")

# Thin elements
@inline thin_bend_no_field(tm, bunch, bendparams)          = error("Undefined for tracking method $tm")                      
@inline thin_bend_bsolenoid(tm, bunch, bendparams, bm0)    = error("Undefined for tracking method $tm")                    
@inline thin_bend_bdipole(tm, bunch, bendparams, bm1)      = error("Undefined for tracking method $tm")                  
@inline thin_bend_bquadrupole(tm, bunch, bendparams, bm2)  = error("Undefined for tracking method $tm")                      
@inline thin_bend_bmultipole(tm, bunch, bendparams, bdict) = error("Undefined for tracking method $tm")                                           
@inline thin_bsolenoid(tm, bunch, bm0)                     = error("Undefined for tracking method $tm")   
@inline thin_bdipole(tm, bunch, bm1)                       = error("Undefined for tracking method $tm") 
@inline thin_bquadrupole(tm, bunch, bm2)                   = error("Undefined for tracking method $tm")     
@inline thin_bmultipole(tm, bunch, bdict)                  = error("Undefined for tracking method $tm")

# Thick elements
@inline thick_bend_no_field(tm, bunch, bendparams, L)          = error("Undefined for tracking method $tm")                      
@inline thick_bend_bsolenoid(tm, bunch, bendparams, bm0, L)    = error("Undefined for tracking method $tm")                   
@inline thick_bend_bdipole(tm, bunch, bendparams, bm1, L)      = error("Undefined for tracking method $tm")                 
@inline thick_bend_bquadrupole(tm, bunch, bendparams, bm2, L)  = error("Undefined for tracking method $tm")                     
@inline thick_bend_bmultipole(tm, bunch, bendparams, bdict, L) = error("Undefined for tracking method $tm")                                          
@inline thick_bsolenoid(tm, bunch, bm0, L)                     = error("Undefined for tracking method $tm")  
@inline thick_bdipole(tm, bunch, bm1, L)                       = error("Undefined for tracking method $tm")
@inline thick_bquadrupole(tm, bunch, bm2, L)                   = error("Undefined for tracking method $tm")    
@inline thick_bmultipole(tm, bunch, bdict, L)                  = error("Undefined for tracking method $tm") 
@inline drift(tm, bunch, L)                                    = error("Undefined for tracking method $tm") 

# === Determine thick or thin === #
@inline bend_no_field(tm, bunch, bendparams, L)          = L ≈ 0 ? thin_bend_no_field(tm, bunch, bendparams)          : thick_bend_no_field(tm, bunch, bendparams, L)          
@inline bend_bsolenoid(tm, bunch, bendparams, bm0, L)    = L ≈ 0 ? thin_bend_bsolenoid(tm, bunch, bendparams, bm0)    : thick_bend_bsolenoid(tm, bunch, bendparams, bm0, L)                       
@inline bend_bdipole(tm, bunch, bendparams, bm1, L)      = L ≈ 0 ? thin_bend_bdipole(tm, bunch, bendparams, bm1)      : thick_bend_bdipole(tm, bunch, bendparams, bm1, L)                       
@inline bend_bquadrupole(tm, bunch, bendparams, bm2, L)  = L ≈ 0 ? thin_bend_bquadrupole(tm, bunch, bendparams, bm2)  : thick_bend_bquadrupole(tm, bunch, bendparams, bm2, L)                       
@inline bend_bmultipole(tm, bunch, bendparams, bdict, L) = L ≈ 0 ? thin_bend_bmultipole(tm, bunch, bendparams, bdict) : thick_bend_bmultipole(tm, bunch, bendparams, bdict, L)                                           
@inline bsolenoid(tm, bunch, bm0, L)                     = L ≈ 0 ? thin_bsolenoid(tm, bunch, bm0)                     : thick_bsolenoid(tm, bunch, bm0, L)                       
@inline bdipole(tm, bunch, bm1, L)                       = L ≈ 0 ? thin_bdipole(tm, bunch, bm1)                       : thick_bdipole(tm, bunch, bm1, L)                       
@inline bquadrupole(tm, bunch, bm2, L)                   = L ≈ 0 ? thin_bquadrupole(tm, bunch, bm2)                   : thick_bquadrupole(tm, bunch, bm2, L)                       
@inline bmultipole(tm, bunch, bdict, L)                  = L ≈ 0 ? thin_bmultipole(tm, bunch, bdict)                  : thick_bmultipole(tm, bunch, bdict, L)                       

