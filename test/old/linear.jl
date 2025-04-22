
using BeamTracking
using Test
using GTPSA


@testset "Linear" begin
  #elements
  #drift
  drift_ele = Linear.Drift(L=0.50)
  #Quadrupole
  quad_ele1 = Linear.Quadrupole(L=0.4, Bn1 = 0.0020) #positive Bn1, quad defocusing for e-
  quad_ele2 = Linear.Quadrupole(L=0.4, Bn1 = -0.0020) #negative Bn1, quad focusing for e-
  quad_ele3 = Linear.Quadrupole(L=0.4, Bn1 = 2.0) #big Bn1
  quad_ele4 = Linear.Quadrupole(L=0.4, Bn1 = 0.0) #no strength, same as drift

  #Sbend
  #sb_ele1 = Linear.SBend(L = 1.5, B0 = 0.000446) # positive field

  #single particle bunch
  d = Descriptor(6,1)
  bi = Bunch(d)
  bf_d = Bunch(d)
  bf_q1 = Bunch(d)
  bf_q2 = Bunch(d)
  bf_q3 = Bunch(d)
  bf_q4 = Bunch(d)
  bf_s1 = Bunch(d)

  #Drift 
  track!(bf_d, drift_ele, bi)
  @test !(bf_d == bi)
  par_d = Particle(bf_d)
  mat_d = GTPSA.jacobian(par_d.z)
  println("\n","drift")
  show(stdout, "text/plain", mat_d)



  #Quad Defocusing for electron
  track!(bf_q1, quad_ele1, bi) # Positive Bn1
  @test !(bf_q1 == bi)
  par_q1 = Particle(bf_q1)
  mat_q1 = GTPSA.jacobian(par_q1.z)
  println("\n","quad with k1=-1.2, Defocusing e-")
  show(stdout, "text/plain", mat_q1)

  #Quad Focusing 
  track!(bf_q2, quad_ele2, bi) # negative Bn1
  @test !(bf_q2 == bi)
  par_q2 = Particle(bf_q2)
  mat_q2 = GTPSA.jacobian(par_q2.z)
  println("\n","quad with k1=1.2, Focusing e-")
  show(stdout, "text/plain", mat_q2)

  #Large k, defocusing
  track!(bf_q3, quad_ele3, bi) # negative Bn1
  @test !(bf_q3 == bi)
  par_q3 = Particle(bf_q3)
  mat_q3 = GTPSA.jacobian(par_q3.z)
  println("\n","quad with B=2, defocusing e-")
  show(stdout, "text/plain", mat_q3)

  #0 strength quadrupole should give the same result as a drift
  track!(bf_q4, quad_ele4, bi) # negative Bn1
  @test !(bf_q4 == bi)
  par_q4 = Particle(bf_q4)
  mat_q4 = GTPSA.jacobian(par_q4.z)
  println("\n","quad with Bn1=0, simplify to drift")
  show(stdout, "text/plain", mat_q4)

  #Sbend
  # track!(sb_ele1, bf_s1, bi) # negative Bn1
  #  @test !(bf_s1 == bi)
  #  par_s1 = Particle(bf_s1)
  #  mat_s1 = GTPSA.jacobian(par_s1.z)
  #  println("\n","Sbend,positive B0")
  #  show(stdout, "text/plain", mat_s1)

end