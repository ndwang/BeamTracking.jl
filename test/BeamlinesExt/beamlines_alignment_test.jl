include("../lattices/alignment_lat.jl")

v1 = [0.1  0.2  0.3  0.4  0.5  0.6]
drift_out = [0.3603778219616477 0.20000000000000004 0.8207556439232953 0.4 0.5565025848905005 0.6]
bend1_out = [0.46408077686019533 0.25826924569290155 0.9250980107869189 0.46372657485893526 0.39699731659260634 0.6]
kicker1_out = [0.35583409051960385 0.19007681187771291 0.8203003409845823 0.399362770761857 0.5571770838946251 0.6]

@testset "BeamlinesAlignment" begin
  b1 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b1, bline_drift1)
  @test b1.coords.v ≈ drift_out

  b2 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b2, bline_drift2)
  @test b2.coords.v ≈ drift_out

  b3 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b3, bline_bend1)
  @test b3.coords.v ≈ bend1_out

  b4 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b4, bline_bend2)
  @test b4.coords.v ≈ drift_out

  b5 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b5, bline_kicker1)
  @test b5.coords.v ≈ kicker1_out
end