#

@inline function aperture(tm, bunch, alignmentP, entering::Bool)
  if !isactive(alignmentP); return nothing; end

  x1 = alignmentP.x1_limit
  x2 = alignmentP.x2_limit
  y1 = alignmentP.y1_limit
  y2 = alignmentP.y2_limit

  if alignmentP.aperture_shape == ApertureShape.Elliptical && (x1 == -Inf || x2 == Inf || y1 == -Inf || y2 == Inf)
    error("Not all x1_limit, x2_limit, y1_limit, y2_limit parameters for elliptical aperture set.")
  end

  if entering
    if alignmentP.aperture_at != ApertureAt.Entrance && alignmentP.aperture_at != ApertureAt.BothEnds
      return nothing
    end
  else
    if alignmentP.aperture_at != ApertureAt.Exit && alignmentP.aperture_at != ApertureAt.BothEnds
      return nothing
    end
  end

  #

  if alignmentP.aperture_shape == ApertureShape.Rectangular
    return KernelCall(BeamTracking.track_aperture_rectangular!, (x1, x2, y1, y2))
  else
    return KernelCall(BeamTracking.track_aperture_elliptical!, (x1, x2, y1, y2))
  end
end
