function track_sigma(v0, line; N=10^6)
  zs = []
  for _ in 1:N 
    b0 = Bunch(deepcopy(v0), R_ref=line.R_ref, species=line.species_ref)
    track!(b0,line)
    push!(zs, [b0.coords.v[i] for i in 1:6])
  end

  centroid = [sum([z[i] for z in zs])/length(zs) for i in 1:6]

  sigma = zeros(6,6)
  for i in 1:6
    for j in 1:6
      sigma[i,j] = sum([(z[i]-centroid[i]) .* (z[j]-centroid[j]) for z in zs])/length(zs)
    end
  end

  return centroid, sigma 
end

bend = SBend(g = 0.01, L = 2.0, 
tracking_method = SplitIntegration(order = 2, num_steps = 10, 
radiation_damping_on = true, radiation_fluctuations_on = true))
line = Beamline([bend], species_ref = Species("electron"), E_ref = 18.0e9)

centroid, sigma = track_sigma(zeros(6), line)