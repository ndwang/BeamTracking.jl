using Test
using GTPSA, BeamTracking
using BeamTracking
using Beamlines
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI, BunchView

function track!(bunch::Bunch, ele::LineElement)
    b = BeamTracking.BunchView(bunch)
    theta = ele.g * ele.L
    mc2 = BeamTracking.massof(bunch.species)
    p0c = 10E6
    tilde_m = mc2/p0c
    beta_0 = p0c / sqrt(mc2^2 + p0c^2)
    return BeamTracking.ExactTracking.exact_bend!(1, b::BunchView, theta, ele.g, ele.K0, tilde_m, beta_0, ele.L)
end

xi  = 0.1
pxi = -7.5e-4
yi  = 0.1
pyi = -3e-4 	
zi  = 0.1
pzi = -1e-3
ps = [xi;; pxi;; yi;; pyi;; zi;; pzi]

d = Descriptor(6, 4)
Δx = transpose(@vars(d))
ps += Δx

Brho = 10E6/BeamTracking.C_LIGHT
b = BeamTracking.Bunch(ps, species = ELECTRON, Brho_ref = Brho)
ele = SBend(L = 0.5, g = 1.0, e1 = 0.0, e2 = 0.0, K0 = 1.001, tracking_method = Exact())
track!(b, ele)
M = GTPSA.jacobian(b.v)

bmadstd = [  0.85374615     0.52189495     0.00000000     0.00004425     0.00000000     0.14734148
            -0.47990496     0.87794249     0.00000000     0.00014397     0.00000000     0.47942570   
            -0.00014415    -0.00003629     1.00000000     0.54837870     0.00000000     0.00015870
             0.00000000     0.00000000     0.00000000     1.00000000     0.00000000     0.00000000   
            -0.48001775    -0.12085250     0.00000000     0.00015870     1.00000000    -0.01861516 
             0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.00000000]

@test M ≈ bmadstd