#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix
# (equal to number of temporaries needed for a single particle)
struct Exact end

MAX_TEMPS(::Exact) = 5

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI
const TRACKING_METHOD = Exact

# Update the reference energy of the canonical coordinates
# BUG: z and pz are not updated correctly
#=
@inline function update_P0!(i, v, work, Brho_initial, Brho_final)
  @inbounds begin
    @FastGTPSA! v[i,PXI] = v[i,PXI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PYI] = v[i,PYI] * Brho_initial / Brho_final
    @FastGTPSA! v[i,PZI] = v[i,PZI] * Brho_initial / Brho_final
  end
  return v
end
=#

# Misalignments (TO-DO: rotational misalignments)
@inline function misalign!(i, v, work, x_offset, y_offset, sgn) #x_rot, y_rot, tilt,
  @assert sgn == -1 || sgn == 1 "Incorrect value for sgn (use -1 if entering, 1 if exiting)"
  @inbounds begin
    @FastGTPSA! v[i,XI] += sgn*x_offset
    @FastGTPSA! v[i,YI] += sgn*y_offset
  end
  return v
end


#
# ===============  E X A C T   D R I F T  ===============
#
"""
exact_drift!()

In the computation of z_final, we use the fact that
    1/√a - 1/√b == (b - a)/(√a √b (√a + √b))
to avoid the potential for severe cancellation when
a and b both have the form 1 + ε for different small
values of ε.
"""
@inline function exact_drift!(i, v, work, tilde_m, gamsqr_0, beta_0, L)
  @assert size(work, 2) >= 1 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 1) for exact_drift!()."
  @inbounds begin @FastGTPSA! begin
    work[i,1] = sqrt((1.0 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
    v[i,XI]   = v[i,XI] + v[i,PXI] * L / work[i,1]
    v[i,YI]   = v[i,YI] + v[i,PYI] * L / work[i,1]
    # high-precision computation of z_final
    # vf.z = vi.z - (1 + δ) * L * (1 / Ps - 1 / (β0 * sqrt((1 + δ)^2 + tilde_m^2)))
    v[i,ZI]   = v[i,ZI] - ( (1.0 + v[i,PZI]) * L
                  * ((v[i,PXI]^2 + v[i,PYI]^2) - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                  / ( beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) * work[i,1]
                      * (beta_0 * sqrt((1 + v[i,PZI])^2 + tilde_m^2) + work[i,1])
                    )
                )
  end end
  return v
end # function exact_drift!()


#
# ===============  Q U A D R U P O L E  ===============
#
"""
mkm_quadrupole!()

This integrator uses the so-called Matrix-Kick-Matrix method to implement
a quadrupole integrator accurate though second-order in the integration
step-size.

beta_gamma_0: reference value of βγ
k2_num: g / Bρ0 = g / (p0 / q)
        where g and Bρ0 respectively denote the quadrupole gradient and
        (signed) reference magnetic rigidity.
"""
@inline function mkm_quadrupole!(i, v, work, beta_gamma_0, k2_num, L)
  @assert size(work, 2) >= 5 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 5) for mkm_quadrupole!()."
  @inbounds begin @FastGTPSA! begin
    quadrupole_matrix!(i, v, work, k2_num, L / 2)
    quadrupole_kick!(i, v, work, beta_gamma_0, L)
    quadrupole_matrix!(i, v, work, k2_num, L / 2)
  end end
  return v
end # function mkm_quadrupole!()


"""
quadrupole_matrix!()

track "matrix part" of a quadrupole
"""
@inline function quadrupole_matrix!(i, v, work, k2_num, s)
  @assert size(work, 2) >= 5 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 5) for quadrupole_matrix!()."
  @inbounds begin @FastGTPSA! begin
    work[i,1] = v[i,PXI] / (1.0 + v[i,PZI])  # x'
    work[i,2] = v[i,PYI] / (1.0 + v[i,PZI])  # y'
    # the following line requires complex arithmetic
    work[i,3] = sqrt(k2_num / (1.0 + v[i,PZI])) * s  # κs for each particle
    work[i,4] = v[i,XI]  * cos(work[i,3])  +        s * work[i,1] * sincu(work[i,3])
    v[i,PXI]  = v[i,PXI] * cos(work[i,3])  - k2_num * s * v[i,XI] * sincu(work[i,3])
    work[i,5] = v[i,YI]  * cosh(work[i,3]) +        s * work[i,2] * sinhcu(work[i,3])
    v[i,PYI]  = v[i,PYI] * cosh(work[i,3]) + k2_num * s * v[i,YI] * sinhcu(work[i,3])
    v[i,ZI]   = (v[i,ZI] - (s / 4) * ( work[i,1]^2 * (1 + sincu(2.0work[i,3])) + work[i,2]^2 * (1 + sinhcu(2.0work[i,3]))
                                       + k2_num / (1.0 + v[i,PZI])
                                         * (v[i,XI]^2 * (1 - sincu(2.0work[i,3])) - v[i,YI]^2 * (1 - sinhcu(2.0work[i,3]))) )
                         + (v[i,XI] * work[i,1] * sin(work[i,3])^2 - v[i,YI] * work[i,2] * sinh(work[i,3])^2) / 2)
    v[i,XI] = work[i,4]
    v[i,YI] = work[i,5]
  end end
  return v
end # function quadrupole_matrix!()


"""
quadrupole_kick!()

track "remaining part" of quadrupole, a position kick

### Implementation
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.
"""
@inline function quadrupole_kick!(i, v, work, beta_gamma_0, s)
  @assert size(work, 2) >= 3 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 3) for quadrupole_kick!"
  @inbounds begin @FastGTPSA! begin
  tilde_m  = 1 / beta_gamma_0  # mc^2 / p0·c
  gamsqr_0 = 1 + beta_gamma_0^2
  beta_0   = beta_gamma_0 / sqrt(gamsqr_0)
  work[i,1] = 1 + v[i,PZI]             # reduced momentum, P/P0 = 1 + δ
  work[i,2] = v[i,PXI]^2 + v[i,PYI]^2  # P⟂^2
  work[i,3] = sqrt(p^2 - ptr2)         # Ps
  v[i,XI] = v[i,XI] + s * v[i,PXI] / work[i,1] * work[i,2] / (work[i,3] * (work[i,1] + work[i,3]))
  v[i,YI] = v[i,YI] + s * v[i,PYI] / work[i,1] * work[i,2] / (work[i,3] * (work[i,1] + work[i,3]))
  v[i,ZI] = v[i,ZI] - s * ( (1.0 + v[i,PZI])
                              * (work[i,2] - v[i,PZI] * (2 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt((1.0 + v[i,PZI])^2 + tilde_m^2) * work[i,3]
                                    * (beta_0 * sqrt((1.0 + v[i,PZI])^2 + tilde_m^2) + work[i,3])
                                  )
                            - work[i,2] / (2 * (1 + v[i,PZI])^2)
                          )
  end end
  return v
end # function quadrupole_kick!()


#
# ===============  M U L T I P O L E  ===============
#
"""
dkd_multipole()

This integrator uses Drift-Kick-Drift to track a beam through
a straight, finite-length multipole magnet. This method is
accurate through second order in the step size. The vectors
bm and am contain the normal and skew multipole strengths, starting
with the dipole component. (For example, b[3] denotes the normal
sextupole strength in Tesla/m^2.) The parameter ns denotes the

Arguments
---------
beta_gamma_0: reference value of βγ
kn: normal multipole strengths
ks: skew multipole strengths
L:  element length
"""
@inline function dkd_multipole!(i, v, work, beta_gamma_0, b_rho, kn, ks, L)
  @assert size(work, 2) >= 5 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 5) for dkd_multipole!()."
  @inbounds begin @FastGTPSA! begin
  tilde_m  = 1 / beta_gamma_0  # mc^2 / p0·c
  gamsqr_0 = 1 + beta_gamma_0^2
  beta_0   = beta_gamma_0 / sqrt(gamsqr_0)
  # ds = L / ns
  # for i = 1:ns
  exact_drift!(i, v, work, tilde_m, gamsqr_0, beta_0, L / 2)
  multipole_kick!(i, v, work, mm, kn * L, ks * L)
  exact_drift!(i, v, work, tilde_m, gamsqr_0, beta_0, L / 2)
  # end
  end end
  return v
end # function dkd_multipole!()


@inline function multipole_kick!(i, v, work, mm, knl, ksl)
  """
  This function tracks a beam of particles through a thin-lens
  multipole having integrated normal and skew strengths in the
  coefficient vectors knl and ksl respectively.
  NP: This function IGNORES the m = 1 (dipole) components.
  """
  @assert size(work, 2) >= 2 && size(work, 1) == N_particle "Size of work matrix must be at least ($N_particle, 2) for multipole_kick!()."
  @inbounds begin @FastGTPSA! begin
  [ work[i,1], work[i,2] ] = dpxy_multipole(mm, knl, ksl, v[i,XI], v[i,YI])
  v[i,PXI] = v[i,PXI] + work[i,1]
  v[i,PYI] = v[i,PYI] + work[i,2]
  end end
  return v
end # function multipole_kick!()


function dpxy_multipole(mm, knl, ksl, x, y)
  """
  This function uses a Horner-like scheme (see Shachinger and
  Talman [SSC-52]) to compute the transverse kicks induced by
  a pure multipole magnet. The algorithm takes advantage of
  the complex representation of the vector potential Az:
   - Re{ Σ_m (bm + i am) (x + i y)^m / m }.
  This method supposedly has good numerical properties, though
  I've not seen a proof of that.

  NB: This function IGNORES the m = 1 (dipole) components of
  both knl and ksl. Moreover, it should *not* include an m = 2
  (quadrupole) component unless this function is used for a
  thin-lens (zero-length) element.

  Arguments
  ---------
  mm:  maximum order, also length of both knl and ksl having
         non-zero coefficients
  knl: vector of normal integrated multipole strengths
  ksl: vector of skew integrated multipole strengths
       NB: Here the m-th component of knl (ksl) denotes the
         normal (skew) component of the m-th multipole after
         scaling by the reference Bρ. For example, knl[2]
         denotes the normal integrated quadrupole strength
         scaled by Bρo.
  x:   horizontal particle coordinate
  y:   vertical particle coordinate
  """
  ar = knl[m] * x - ksl[m] * y
  ai = knl[m] * y + ksl[m] * x
  while m > 2,
    m -= 1
    t  = (knl[m] * x - ksl[m] * y) + (ar * x - ai * y)
    ai = (knl[m] * y + ksl[m] * x) + (ar * y + ai * x)
    ar = t
  end
  return [ -ar, ai ]
end # function dpxy_multipole()


function binom(m::Integer, x, y)
  """
  This function computes the real and imaginary parts of
  (x + i y)^m. One can use these components to compute,
  among other things, the multipole kick induced by two-
  dimensional multipole magnets.
  """
  if m == 0
    return [ 1.0, 0.0 ]
  end
  ar = x
  ai = y
  mm = m
  while mm > 1
    mm -= 1
    t  = x * ar - y * ai
    ai = y * ar + x * ai
    ar = t
  end
  return [ ar, ai ]
end # function binom()


end
