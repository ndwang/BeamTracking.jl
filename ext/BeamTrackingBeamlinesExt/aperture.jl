#

@inline function aperture(tm, bunch, alignmentP, entering::Bool)
  if !isactive(alignmentP); return nothing; end

  x1 = alignmentP.x1_limit
  x2 = alignmentP.x2_limit
  y1 = alignmentP.y1_limit
  y2 = alignmentP.y2_limit

  if alignmentP.aperture_shape == ApertureShape.elliptical && (x1 == -Inf || x2 == Inf || y1 == -Inf || y2 == Inf)
    return nothing
  end

  if entering
    if alignmentP.aperture_at != ApertureAt.Entrance && alignmentP.aperture_at != BothEnds
      return nothing
    end
  else
    if alignmentP.aperture_at != ApertureAt.Exit && alignmentP.aperture_at != BothEnds
      return nothing
    end
  end

  #

  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)

  if alignmentP.aperture_shape == ApertureShape.Rectangular
    return KernelCall(track_aperture_rectangular!, (x1, x2, y1, y2))
  else
    return KernelCall(track_aperture_elliptical!, (x1, x2, y1, y2))
  end
end
