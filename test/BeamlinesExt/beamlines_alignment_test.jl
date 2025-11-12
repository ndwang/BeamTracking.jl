include("../lattices/alignment_lat.jl")

v1 = [0.1  0.2  0.3  0.4  0.5  0.6]
drift_out = [0.3603778219616477 0.20000000000000004 0.8207556439232953 0.4 0.5565025848905005 0.6]
bend1_out = [0.46408077686019533 0.25826924569290155 0.9250980107869189 0.46372657485893526 0.39699731659260634 0.6]
v3 = [0.0 0.0 0.0 0.0 0.0 0.0]

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
end