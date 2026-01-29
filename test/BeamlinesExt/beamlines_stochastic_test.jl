using Random 

@testset "Stochastic radiation" begin
  Random.seed!(0)

  p_over_q_ref = BeamTracking.E_to_R(Species("electron"), 18e9)
  bend = SBend(g = 0.01, L = 2.0, 
  tracking_method = Yoshida(order = 2, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true))
  line = Beamline([bend], species_ref = Species("electron"), p_over_q_ref = p_over_q_ref)

  v0 = [0.01 0.02 0.03 0.04 0.05 0.06]
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889165969752571 
  0.021174046860678353 
  0.10556186702264685 
  0.03999912146491459 
  0.04761050793826471 
  0.05997673560385752]'

  bend.tracking_method = Yoshida(order = 2, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889193258818915
  0.021174302497216673
  0.10556186827683302
  0.03999908434104034
  0.04761050068397177
  0.059975733285604405]'

  bend.tracking_method = Yoshida(order = 4, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.048891974265822244
  0.02117431265208611
  0.10556186796033079
  0.03999910506374448
  0.047610499467786394
  0.059976284409328576]'

  bend.tracking_method = Yoshida(order = 4, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.0488920016565342
  0.021174218561770715
  0.10556186797663822
  0.039998878356961864
  0.04761049873848289
  0.05997027467880636]'


  bend.tracking_method = Yoshida(order = 6, num_steps = 1,
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈ 
  [0.04889226281541181
  0.021174668413673003
  0.10556186827012498
  0.03999920072349538
  0.04761049187637356
  0.059978810393760046]'


  bend.tracking_method = Yoshida(order = 6, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)
  
  @test b0.coords.v ≈
  [0.04889189599231092
  0.02117408337833361
  0.10556186775984797
  0.03999890918412512
  0.04761050140381129
  0.059971097099417885]'


  bend.tracking_method = Yoshida(order = 8, num_steps = 1, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈
  [0.04889209188518841
  0.02117482934639319
  0.10556186796145814
  0.039999846102478705
  0.047610496377937794
  0.059995919985280796]'


  bend.tracking_method = Yoshida(order = 8, num_steps = 2, 
  radiation_damping_on = true, radiation_fluctuations_on = true)
  b0 = Bunch(copy(v0), species = line.species_ref, p_over_q_ref = line.p_over_q_ref)
  track!(b0, line)

  @test b0.coords.v ≈
  [0.04889219169051286
  0.02117499046271834
  0.10556186830041747
  0.039999857154632036
  0.04761049389108041
  0.05999621385138163]'

end
 
