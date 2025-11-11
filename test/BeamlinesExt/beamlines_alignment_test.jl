include("../lattices/alignment_lat.jl")

v1 = [0.1  0.2  0.3  0.4  0.5  0.6]
v2 = [0.3603778219616477 0.20000000000000004 0.8207556439232953 0.4 -3.5830243461386777 0.6]

@testset "BeamlinesAlignment" begin
  b1 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b1, bline_d1)
  @test b1.coords.v ≈ v2

  b2 = Bunch(deepcopy(v1), species=Species("electron"))
  track!(b2, bline_d1)
  @test b2.coords.v ≈ v2
end