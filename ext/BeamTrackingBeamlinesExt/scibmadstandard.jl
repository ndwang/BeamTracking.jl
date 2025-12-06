_track!(
  coords::Coords,
  bunch::Bunch,
  ele::LineElement, 
  tm::SciBmadStandard,
  ramp_without_rf;
  kwargs...
) = _track!(coords, bunch, ele, SplitIntegration(), ramp_without_rf; kwargs...)