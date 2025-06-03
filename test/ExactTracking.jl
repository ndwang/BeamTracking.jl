max_temps = BeamTracking.MAX_TEMPS(Exact())
c_light = 299792458.0
degree = π / 180

# define element paameterss
# -- drifts
ld1 = 0.15;  # m
ld2 = 0.75;  # m
ld3 = 2.00;  # m
ld4 = 2.00;  # m

# -- quadrupoles
lq1 = 0.05;  # m
lq2 = 0.30;  # m
lq3 = 0.75;  # m
lq4 = 1.20;  # m
gr1 = 0.20;  # T/m
gr2 = 0.60;  # T/m
gr3 = 1.20;  # T/m
gr4 = 3.50;  # T/m

# -- multipoles
lm1 = 0.05;  # m
lm2 = 0.30;  # m
lm3 = 0.75;  # m
lm4 = 1.20;  # m
mm1 = [ 3, 4, 5, 6 ]
bv1 = [ 0.4, 0.3, 0.2, 0.1 ];  # T·m^{n-1}

# define beams
# -- species
#e_minus = Species("electron")
#p_plus =  Species("proton")
#mec2 = massof(e_minus) # 0.51099895069 MeV
#mpc2 = massof(p_plus)  # 938.27208943 MeV
mec2 = 0.51099895069e6 # eV
mpc2 = 938.27208943e6  # eV
# -- kinetic energy
ek1 =   5.e3;  # eV
ek2 =   1.e6;  # eV
ek3 =   1.e9;  # eV
ek4 = 250.e9;  # eV
# -- βγ
βγ1 = sqrt(ek1 / mec2 * (ek1 / mec2 + 2))
βγ2 = sqrt(ek2 / mec2 * (ek2 / mec2 + 2))
βγ3 = sqrt(ek3 / mec2 * (ek3 / mec2 + 2))
βγ4 = sqrt(ek4 / mpc2 * (ek4 / mpc2 + 2))
# -- γ^2
γsq1 = 1.0 + βγ1^2
γsq2 = 1.0 + βγ2^2
γsq3 = 1.0 + βγ3^2
γsq4 = 1.0 + βγ4^2
# -- β
β1 = βγ1 / sqrt(γsq1)
β2 = βγ2 / sqrt(γsq2)
β3 = βγ3 / sqrt(γsq3)
β4 = βγ4 / sqrt(γsq4)
# -- pc
pc1 = mec2 * βγ1
pc2 = mec2 * βγ2
pc3 = mec2 * βγ3
pc4 = mpc2 * βγ4
# -- Bρ
Bρ1 = -pc1 / c_light
Bρ2 = -pc2 / c_light
Bρ3 = -pc3 / c_light
Bρ4 =  pc4 / c_light

# -- beams
# beam1.initial
xi  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ]
pxi = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ]
yi  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ]
pyi = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ]
zi  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ]
pzi = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ]

# beam1.dr1.final
xf_dr1  = [ 0.,  0.,                     0.,                     2.e-3,                  2.1123876489808660e-3, -2.1123876489808660e-3 ]
pxf_dr1 = [ 0.,  0.,                     0.,                     0.,                     7.5e-4,                -7.5e-4                ]
yf_dr1  = [ 0.,  0.,                     0.,                     1.e-3,                  1.0449550595923462e-3, -1.0449550595923462e-3 ]
pyf_dr1 = [ 0.,  0.,                     0.,                     0.,                     3.e-4,                 -3.e-4                 ]
zf_dr1  = [ 0.,  1.4710284465059844e-4, -1.4711135596840458e-4,  1.4710284465059844e-4,  1.4705400485512816e-4,  1.4705400485512816e-4 ]
pzf_dr1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.dr2.final
xf_dr2  = [ 0., 0., 0., 0.002, 0.0025619382449043287, -0.0025619382449043287 ]
pxf_dr2 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr2  = [ 0., 0., 0., 0.001, 0.0012247752979617315, -0.0012247752979617315 ]
pyf_dr2 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr2  = [ 0., 0.00008566359457101641, -0.00008589149602558208, 0.00008566359457101641, 0.00008541939559366522, 0.00008541939559366522 ]
pzf_dr2 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam3.dr3.final
xf_dr3  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
pxf_dr3 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr3  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
pyf_dr3 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr3  = [ 0., 5.209250185095532e-10, -5.224901403178192e-10, 5.209250185095532e-10, -6.506763479180273e-7, -6.506763479180273e-7 ]
pzf_dr3 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam4.dr4.final
xf_dr4  = [ 0., 0., 0., 0.002, 0.003498501986411544,  -0.003498501986411544 ]
pxf_dr4 = [ 0., 0., 0., 0.,    0.00075,               -0.00075 ]
yf_dr4  = [ 0., 0., 0., 0.001, 0.0015994007945646172, -0.0015994007945646172 ]
pyf_dr4 = [ 0., 0., 0., 0.,    0.0003,                -0.0003 ]
zf_dr4  = [ 0., 2.7919184691863886e-8, -2.800306686850912e-8, 2.7919184691863886e-8, -6.232780882446728e-7, -6.232780882446728e-7 ]
pzf_dr4 = [ 0., 0.001, -0.001, 0.001, 0.001, 0.001 ]

# beam1.qf1.final
xf_qf1  = [ 0.,  0.,                      0.,                     4.4834792779340600e-3,  4.5356143287504990e-3, -4.5356143287504990e-3 ]
pxf_qf1 = [ 0.,  0.,                      0.,                     1.1607821136240948e-1,  1.1776162121669208e-1, -1.1776162121669208e-1 ]
yf_qf1  = [ 0.,  0.,                      0.,                     1.2400905948673489e-4,  1.3427609296030678e-4, -1.3427609296030678e-4 ]
pyf_qf1 = [ 0.,  0.,                      0.,                    -2.8691666098954356e-2, -2.8653744321335432e-2,  2.8653744321335432e-2 ]
zf_qf1  = [ 0.,  4.903428155019947e-5,   -4.903711865613486e-5,  -4.8701323139842656e-5, -5.1605970700562340e-5, -5.1605970700562340e-5 ]
pzf_qf1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qd1.final
xf_qd1  = [ 0.,  0.,                      0.,                     2.4834727340278294e-4,  2.7409827330240064e-4, -2.7409827330240064e-4 ]
pxf_qd1 = [ 0.,  0.,                      0.,                    -5.7391734255615580e-2, -5.7299059109537503e-2,  5.7299059109537503e-2 ]
yf_qd1  = [ 0.,  0.,                      0.,                     2.2414071638535435e-3,  2.2621899982980444e-3, -2.2621899982980444e-3 ]
pyf_qd1 = [ 0.,  0.,                      0.,                     5.8033153149417160e-2,  5.8705242601074584e-2, -5.8705242601074584e-2 ]
zf_qd1  = [ 0.,  4.9034281550199470e-5,  -4.9037118656134860e-5, -1.1242623339116514e-5, -1.1118475705649755e-5, -1.1118475705649755e-5 ]
pzf_qd1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.qf2.final
xf_qf2 =  [ 0.,  0.,                      0.,                     2.9271671401041872e-2,  3.0251479090608963e-2, -3.0251479090608963e-2 ]
pxf_qf2 = [ 0.,  0.,                      0.,                     3.2854870071095094e-1,  3.3959282727762820e-1, -3.3959282727762820e-1 ]
yf_qf2 =  [ 0.,  0.,                      0.,                    -9.7278287676729380e-4, -9.7883210731348450e-4,  9.7883210731348450e-4 ]
pyf_qf2 = [ 0.,  0.,                      0.,                     2.6415505265087600e-3,  2.3544440936341645e-3, -2.3544440936341645e-3 ]
zf_qf2 =  [ 0.,  3.4265437828406564e-5,  -3.4356598410232830e-5, -2.3378443788359830e-3, -2.5012489357825860e-3, -2.5012489357825860e-3 ]
pzf_qf2 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam2.qd2.final
xf_qd2 =  [ 0.,  0.,                      0.,                    -0.0019464186216836795, -0.001961645985092862, 0.001961645985092862 ]
pxf_qd2 = [ 0.,  0.,                      0.,                     0.005200326165715375,   0.004472438532736971, -0.004472438532736971 ]
yf_qd2 =  [ 0.,  0.,                      0.,                     0.014608704236165976,   0.014997970393040412, -0.014997970393040412 ]
pyf_qd2 = [ 0.,  0.,                      0.,                     0.16398929979812307,    0.16837903611031713, -0.16837903611031713 ]
zf_qd2 =  [ 0.,  3.4265437828406564e-5,  -3.435659841023283e-5,  -5.899497183747271e-4,  -6.222650472445767e-4, -6.222650472445767e-4 ]
pzf_qd2 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam3.qf3.final
xf_qf3 =  [ 0.,  0.,                      0.,                     0.002205479704035602,   0.0027865339903883585, -0.0027865339903883585 ]
pxf_qf3 = [ 0.,  0.,                      0.,                     0.0005576983180317967,  0.0013847532608658173, -0.0013847532608658173 ]
yf_qf3 =  [ 0.,  0.,                      0.,                     0.0009006624003320741,  0.0011179443238845692, -0.0011179443238845692 ]
pyf_qf3 = [ 0.,  0.,                      0.,                    -0.0002606852268330615,  9.513485195908057e-6, -9.513485195908057e-6 ]
zf_qf3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.630230762164262e-8,  -4.3665730716563075e-7, -4.3665730716563075e-7 ]
pzf_qf3 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam3.qd3.final
xf_qd3 =  [ 0.,  0.,                      0.,                     0.0018013248008431747,   0.0023445295165179687, -0.0023445295165179687 ]
pxf_qd3 = [ 0.,  0.,                      0.,                    -0.0005213704536906773,   0.0001541263391654702, -0.0001541263391654702 ]
yf_qd3 =  [ 0.,  0.,                      0.,                     0.0011027398519220506,   0.0013351614586315588, -0.0013351614586315588 ]
pyf_qd3 = [ 0.,  0.,                      0.,                     0.00027884915900320065,  0.0006096711218370137, -0.0006096711218370137 ]
zf_qd3 =  [ 0.,  1.9534688194108246e-10, -1.959338026191822e-10, -4.4102076263621526e-8,  -1.6795543457606254e-7, -1.6795543457606254e-7 ]
pzf_qd3 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# beam4.qf4.final
xf_qf4 =  [ 0.,  0.,                      0.,                     0.0019939877700227713,   0.002892187842249561, -0.002892187842249561 ]
pxf_qf4 = [ 0.,  0.,                      0.,                    -0.000010025375230046909, 0.0007377200378071899, -0.0007377200378071899 ]
yf_qf4 =  [ 0.,  0.,                      0.,                     0.0010030091302526826,   0.0013630102695159176, -0.0013630102695159176 ]
pyf_qf4 = [ 0.,  0.,                      0.,                     5.022748546347732e-6,    0.0003059254879156335, -0.0003059254879156335 ]
zf_qf4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726401735327598e-8,  -3.698308917351584e-7, -3.698308917351584e-7 ]
pzf_qf4 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# beam4.qd4.final
xf_qd4 =  [ 0.,  0.,                      0.,                     0.002006018260505365,    0.00290602111426842, -0.00290602111426842 ]
pxf_qd4 = [ 0.,  0.,                      0.,                     0.000010045497092695465, 0.0007623023455299651, -0.0007623023455299651 ]
yf_qd4 =  [ 0.,  0.,                      0.,                     0.0009969938850113856,   0.0013562739161008426, -0.0013562739161008426 ]
pyf_qd4 = [ 0.,  0.,                      0.,                    -5.012687615023454e-6,    0.0002940854775943521, -0.0002940854775943521 ]
zf_qd4 =  [ 0.,  1.675151081511833e-8,   -1.680184012110547e-8,   1.6726365460219405e-8,  -3.78176641902771e-7, -3.78176641902771e-7 ]
pzf_qd4 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                   1.e-3,                  1.e-3                 ]

# beam1.mp1.final
xf_mp1  = [ 0.,  0.,                      0.,                     2.0038906474613948e-3,  2.0414381618753773e-3, -2.0334840740152908e-3 ]
pxf_mp1 = [ 0.,  0.,                      0.,                     1.5578151070769097e-4,  9.0918337810996800e-4, -5.9070183442735890e-4 ]
yf_mp1  = [ 0.,  0.,                      0.,                     9.9028590464704330e-4,  1.0051011842090445e-3, -1.0248329717192377e-3 ]
pyf_mp1 = [ 0.,  0.,                      0.,                    -3.8895234385991640e-4, -9.5748642016826030e-5, -6.9431180204891080e-4 ]
zf_mp1  = [ 0.,  4.9034281550199474e-5,  -4.9037118656134856e-5,  4.903209153457148e-5,   4.901571516330201e-5, 4.9015774852461346e-5   ]
pzf_mp1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.mm1.final
xf_mm1  = [ 0.,  0.,                      0.,                     1.9961093525386053e-3,  2.0334869408519317e-3, -2.041441028695591e-3  ]
pxf_mm1 = [ 0.,  0.,                      0.,                    -1.5578151070769097e-4,  5.9081662189003210e-4, -9.092981655726411e-4  ]
yf_mm1  = [ 0.,  0.,                      0.,                     1.0097140953529567e-3,  1.0248688568666960e-3, -1.005137069338170e-3  ]
pyf_mm1 = [ 0.,  0.,                      0.,                     3.8895234385991640e-4,  6.9574864201682590e-4,  9.431180204891087e-5  ]
zf_mm1  = [ 0.,  4.9034281550199474e-5,  -4.9037118656134856e-5,  4.903209153457148e-5,   4.901574824419741e-5,   4.901571596604123e-5  ]
pzf_mm1 = [ 0.,  1.e-3,                  -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]


# test individual elements
@testset "ExactTracking" begin

  # === drifts ===
  #
  # 5 keV electron
  v = [ xi pxi yi pyi zi pzi ]
  work = zeros(size(v, 1), max_temps)
  BeamTracking.launch!(exact_drift!, v, work, β1, γsq1, 1/βγ1, ld1)
  @test v[:,BeamTracking.XI]  ≈  xf_dr1 (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_dr1 (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_dr1 (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] == pxf_dr1
  @test v[:,BeamTracking.PYI] == pyf_dr1
  @test v[:,BeamTracking.PZI] == pzf_dr1
  #
  # 1 MeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(exact_drift!, v, work, β2, γsq2, 1/βγ2, ld2)
  @test v[:,BeamTracking.XI]  ≈  xf_dr2 (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_dr2 (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_dr2 (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] == pxf_dr2
  @test v[:,BeamTracking.PYI] == pyf_dr2
  @test v[:,BeamTracking.PZI] == pzf_dr2
  #
  # 1 GeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(exact_drift!, v, work, β3, γsq3, 1/βγ3, ld3)
  @test v[:,BeamTracking.XI]  ≈  xf_dr3 (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_dr3 (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_dr3 (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] == pxf_dr3
  @test v[:,BeamTracking.PYI] == pyf_dr3
  @test v[:,BeamTracking.PZI] == pzf_dr3
  #
  # 250 GeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(exact_drift!, v, work, β4, γsq4, 1/βγ4, ld4)
  @test v[:,BeamTracking.XI]  ≈  xf_dr4 (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_dr4 (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_dr4 (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] == pxf_dr4
  @test v[:,BeamTracking.PYI] == pyf_dr4
  @test v[:,BeamTracking.PZI] == pzf_dr4

  # === quadrupoles ===
  #
  # 5 keV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β1, γsq1, 1/βγ1,  gr1 / Bρ1, lq1)
  @test v[:,BeamTracking.XI]  ≈  xf_qf1  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qf1  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qf1  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qf1 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qf1 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qf1
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β1, γsq1, 1/βγ1, -gr1 / Bρ1, lq1)
  @test v[:,BeamTracking.XI]  ≈  xf_qd1  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qd1  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qd1  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qd1 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qd1 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qd1
  #
  # 1 MeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β2, γsq2, 1/βγ2,  gr2 / Bρ2, lq2)
  @test v[:,BeamTracking.XI]  ≈  xf_qf2  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qf2  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qf2  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qf2 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qf2 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qf2
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β2, γsq2, 1/βγ2, -gr2 / Bρ2, lq2)
  @test v[:,BeamTracking.XI]  ≈  xf_qd2  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qd2  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qd2  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qd2 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qd2 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qd2
  #
  # 1 GeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β3, γsq3, 1/βγ3,  gr3 / Bρ3, lq3)
  @test v[:,BeamTracking.XI]  ≈  xf_qf3  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qf3  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qf3  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qf3 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qf3 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qf3
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β3, γsq3, 1/βγ3, -gr3 / Bρ3, lq3)
  @test v[:,BeamTracking.XI]  ≈  xf_qd3  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qd3  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qd3  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qd3 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qd3 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qd3
  #
  # 250 GeV electron
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β4, γsq4, 1/βγ4,  gr4 / Bρ4, lq4)
  @test v[:,BeamTracking.XI]  ≈  xf_qf4  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qf4  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qf4  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qf4 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qf4 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qf4
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(mkm_quadrupole!, v, work, β4, γsq4, 1/βγ4, -gr4 / Bρ4, lq4)
  @test v[:,BeamTracking.XI]  ≈  xf_qd4  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_qd4  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_qd4  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_qd4 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_qd4 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_qd4

#  # === multipoles ===
#  #
#  # 5 keV electron
  v = [ xi pxi yi pyi zi pzi ]
  kn1 = bv1 * cos(15degree) / Bρ1
  ks1 = bv1 * sin(15degree) / Bρ1
  BeamTracking.launch!(dkd_multipole!, v, work, β1, γsq1, 1/βγ1, mm1, kn1, ks1, lm1)
  @test v[:,BeamTracking.XI]  ≈  xf_mp1  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_mp1  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_mp1  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_mp1 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_mp1 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_mp1
  v = [ xi pxi yi pyi zi pzi ]
  BeamTracking.launch!(dkd_multipole!, v, work, β1, γsq1, 1/βγ1, mm1, -kn1, -ks1, lm1)
  @test v[:,BeamTracking.XI]  ≈  xf_mm1  (rtol=5.e-13)
  @test v[:,BeamTracking.YI]  ≈  yf_mm1  (rtol=5.e-13)
  @test v[:,BeamTracking.ZI]  ≈  zf_mm1  (rtol=5.e-13)
  @test v[:,BeamTracking.PXI] ≈  pxf_mm1 (rtol=5.e-13)
  @test v[:,BeamTracking.PYI] ≈  pyf_mm1 (rtol=5.e-13)
  @test v[:,BeamTracking.PZI] == pzf_mm1


end

