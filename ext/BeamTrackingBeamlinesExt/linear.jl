
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]


# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  v,
  work,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::Linear
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L

  # Function barrier
  linear_universal!(i, v, work, bunch, L, bm, bp, ma)
end


# Step 2: Push particles through -----------------------------------------
function linear_universal!(
  i, 
  v, 
  work,
  bunch,
  L, 
  bmultipoleparams, 
  bendparams, 
  alignmentparams
) 
  if isactive(bendparams)
    error("bend tracking not implemented yet")
  end

  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)

  if !isactive(bmultipoleparams)
    runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2)
  else 
    if any(t -> t > 2, keys(bmultipoleparams.bdict)) || !haskey(bmultipoleparams.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end

    bm1 = bmultipoleparams.bdict[2]

    s1 = bm1.strength

    if !bm1.normalized
      s1 /= bunch.Brho_ref
    end

    if L != 0
      if bm1.integrated
        K1 = s1/L
      end
      mx, my = LinearTracking.linear_quad_matrices(K1, L)
    else
      if !bm1.integrated
        error("LineElement length is zero but has a non-integrated magnetic multipole")
      end
      mx, my = LinearTracking.linear_thin_quad_matrices(s1)
    end

    r56 = L/gamma_0^2
    runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56)
  end
  return v
end
