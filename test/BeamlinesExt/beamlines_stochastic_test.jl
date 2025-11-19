using Random 

@testset "Stochastic radiation" begin
  Random.seed!(0)

  R_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  bend = SBend(g = 0.01, L = 2.0, 
  tracking_method = SplitIntegration(order = 2, radiation_damping_on = true, 
  radiation_fluctuations_on = true))
  line = Beamline([bend], species_ref = Species("electron"), R_ref = R_ref)

  v0 = [0.01 0.02 0.03 0.04 0.05 0.06]
  b0 = Bunch(copy(v0))
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.048892166423554095 
  0.02117516850989616 
  0.10556186816301837 
  0.040000228166203126 
  0.0476104945848721 
  0.06000604532567636]'

  bend.tracking_method = SplitIntegration(order = 4, radiation_damping_on = true, 
  radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0))
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889184098942946 
  0.021173991360472644 
  0.10556186772381217 
  0.03999876450274291 
  0.047610502977529556 
  0.05996725973596145]'

  bend.tracking_method = SplitIntegration(order = 6, radiation_damping_on = true, 
  radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0))
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889213616438054 
  0.021174519757519893 
  0.10556186809955524 
  0.03999917292581259 
  0.04761049520962324 
  0.05997807693141073]'

  bend.tracking_method = SplitIntegration(order = 8, radiation_damping_on = true, 
  radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0))
  track!(b0, line)

  @test b0.coords.v ≈
  [0.04889210975355922 
  0.021174436067643392 
  0.10556186804578076 
  0.03999906716683646 
  0.047610495905398675
  0.05997528633046387]'
end
 
