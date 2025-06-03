#=

Exact tracking methods

=#
# Define the Exact tracking method, and number of columns in the work matrix
# (equal to number of temporaries needed for a single particle)
struct Exact end

MAX_TEMPS(::Exact) = 7

module ExactTracking
using ..GTPSA, ..BeamTracking, ..StaticArrays
using ..BeamTracking: XI, PXI, YI, PYI, ZI, PZI
const TRACKING_METHOD = Exact

export exact_drift!
export mkm_quadrupole!, quadrupole_matrix!, quadrupole_kick!
export dkd_multipole!, multipole_kick!
export exact_sbend!

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

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
L: element length
"""
@inline function exact_drift!(i, v, work, beta_0, gamsqr_0, tilde_m, L)
  @assert size(work, 2) >= 1 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 1) for exact_drift!()."
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

This integrator uses Matrix-Kick-Matrix to implement a quadrupole
integrator accurate though second-order in the integration step-size.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
k2_num:   g / Bρ0 = g / (p0 / q)
          where g and Bρ0 respectively denote the quadrupole gradient
          and (signed) reference magnetic rigidity.
L: element length
"""
@inline function mkm_quadrupole!(i, v, work, beta_0, gamsqr_0, tilde_m, k2_num, L)
  @assert size(work, 2) >= 7 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 7) for mkm_quadrupole!()."
  @inbounds begin #@FastGTPSA! begin
    #ds = L / ns
    #for i = 1:ns
    quadrupole_matrix!(i, v, work, k2_num, L / 2)
    quadrupole_kick!(  i, v, work, beta_0, gamsqr_0, tilde_m, L)
    quadrupole_matrix!(i, v, work, k2_num, L / 2)
    #end
  end #end
  return v
end # function mkm_quadrupole!()


"""
quadrupole_matrix!()

Track "matrix part" of quadrupole.

Arguments
—————————
k2_num:  g / Bρ0 = g / (p0 / q)
         where g and Bρ0 respectively denote the quadrupole gradient
         and (signed) reference magnetic rigidity.
s: element length
"""
@inline function quadrupole_matrix!(i, v, work, k2_num, s)
  @assert size(work, 2) >= 7 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 7) for quadrupole_matrix!()."
  @inbounds begin #@FastGTPSA! begin
    sgn = sign(k2_num)
    focus   = k2_num >= 0  # horizontally focusing for positive particles
    defocus = k2_num <  0  # horizontally defocusing for positive particles

    work[i,1] = v[i,PXI] / (1.0 + v[i,PZI])  # x'
    work[i,2] = v[i,PYI] / (1.0 + v[i,PZI])  # y'
    work[i,3] = sqrt(abs(k2_num / (1.0 + v[i,PZI]))) * s  # |κ|s
    work[i,4] = focus * cos(work[i,3])    + defocus * cosh(work[i,3])
    work[i,5] = focus * cosh(work[i,3])   + defocus * cos(work[i,3])
    work[i,6] = focus * sincu(work[i,3])  + defocus * sinhcu(work[i,3])
    work[i,7] = focus * sinhcu(work[i,3]) + defocus * sincu(work[i,3])

    v[i,PXI] = v[i,PXI] * work[i,4] - k2_num * s * v[i,XI] * work[i,6]
    v[i,PYI] = v[i,PYI] * work[i,5] + k2_num * s * v[i,YI] * work[i,7]
    v[i,ZI]  = (v[i,ZI] - (s / 4) * (  work[i,1]^2 * (1.0 + work[i,6] * work[i,4])
                                     + work[i,2]^2 * (1.0 + work[i,7] * work[i,5])
                                     + k2_num / (1.0 + v[i,PZI])
                                         * ( v[i,XI]^2 * (1.0 - work[i,6] * work[i,4])
                                           - v[i,YI]^2 * (1.0 - work[i,7] * work[i,5]) )
                                    )
                        + sgn * ( v[i,XI] * work[i,1] * (work[i,3] * work[i,6])^2
                                - v[i,YI] * work[i,2] * (work[i,3] * work[i,7])^2 ) / 2.0
               )
    v[i,XI]  = v[i,XI] * work[i,4] + work[i,1] * s * work[i,6]
    v[i,YI]  = v[i,YI] * work[i,5] + work[i,2] * s * work[i,7]
  end #end
  return v
end # function quadrupole_matrix!()


"""
quadrupole_kick!()

Track "remaining part" of quadrupole —— a position kick.

### Note re implementation:
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
s: element length
"""
@inline function quadrupole_kick!(i, v, work, beta_0, gamsqr_0, tilde_m, s)
  @assert size(work, 2) >= 3 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 3) for quadrupole_kick!"
  @inbounds begin #@FastGTPSA! begin
  work[i,1] = 1.0 + v[i,PZI]                 # reduced total momentum,  P/P0 = 1 + δ
  work[i,2] = v[i,PXI]^2 + v[i,PYI]^2        # (transverse momentum)^2, P⟂^2 = (Px^2 + Py^2) / P0^2
  work[i,3] = sqrt(work[i,1]^2 - work[i,2])  # longitudinal momentum,   Ps = √[(1 + δ)^2 - P⟂^2]
  v[i,XI] = v[i,XI] + s * v[i,PXI] / work[i,1] * work[i,2] / (work[i,3] * (work[i,1] + work[i,3]))
  v[i,YI] = v[i,YI] + s * v[i,PYI] / work[i,1] * work[i,2] / (work[i,3] * (work[i,1] + work[i,3]))
  v[i,ZI] = v[i,ZI] - s * ( work[i,1]
                              * (work[i,2] - v[i,PZI] * (2.0 + v[i,PZI]) / gamsqr_0)
                                / ( beta_0 * sqrt(work[i,1]^2 + tilde_m^2) * work[i,3]
                                    * (beta_0 * sqrt(work[i,1]^2 + tilde_m^2) + work[i,3])
                                  )
                            - work[i,2] / (2 * work[i,1]^2)
                          )
  end #end
  return v
end # function quadrupole_kick!()


#
# ===============  M U L T I P O L E  ===============
#
"""
dkd_multipole()

This integrator uses Drift-Kick-Drift to track a beam through
a straight, finite-length multipole magnet. The method is
accurate through second order in the step size. The vectors
kn and ks contain the normal and skew multipole strengths,
starting with the dipole component. (For example, b[3] denotes
the normal sextupole strength in Tesla/m^2.) The argument ns
denotes the number of slices.

Arguments
—————————
beta_0:   β_0 = (βγ)_0 / √(γ_0^2)
gamsqr_0: γ_0^2 = 1 + (βγ)_0^2
tilde_m:  1 / (βγ)_0  # mc^2 / p0·c
mm: vector of m values for non-zero multipole coefficients
kn: vector of normal multipole strengths scaled by Bρ0
ks: vector of skew multipole strengths scaled by Bρ0
L:  element length
"""
@inline function dkd_multipole!(i, v, work, beta_0, gamsqr_0, tilde_m, mm, kn, ks, L)
  @assert size(work, 2) >= 3 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 3) for dkd_multipole!()."
  @inbounds begin #@FastGTPSA! begin
    #ds = L / ns
    #for i = 1:ns
    exact_drift!(   i, v, work, beta_0, gamsqr_0, tilde_m, L / 2)
    multipole_kick!(i, v, work, mm, kn * L, ks * L)
    exact_drift!(   i, v, work, beta_0, gamsqr_0, tilde_m, L / 2)
    #end
  end #end
  return v
end # function dkd_multipole!()


"""
    multipole_kick!(i, v, work, mm, knl, ksl)

Track a beam of particles through a thin-lens multipole
having integrated normal and skew strengths listed in the
coefficient vectors knl and ksl respectively. The vector mm
lists the order of the corresponding entries in knl and ksl.

NB: This function IGNORES the m = 1 (dipole) components of
both knl and ksl. Moreover, it should *not* include an m = 2
(quadrupole) component unless this function is used for a
thin-lens (zero-length) element.

The algorithm used in this function takes advantage of the
complex representation of the vector potential Az,
  - ``-Re{ sum_m (b_m + i a_m) (x + i y)^m / m }``,
and uses a Horner-like scheme (see Shachinger and Talman
[SSC-52]) to compute the transverse kicks induced by a pure
multipole magnet.  This method supposedly has good numerical
properties, though I've not seen a proof of that.

DTA: Ordering matters!
DTA: Add thin dipole kick.

### Arguments
 - mm:  vector of m values for non-zero multipole coefficients </br>
 - knl: vector of normal integrated multipole strengths </br>
 - ksl: vector of skew integrated multipole strengths </br>


     NB: Here the j-th component of knl (ksl) denotes the
       normal (skew) component of the multipole strength of
       order mm[j] (after scaling by the reference Bρ).
       For example, if mm[j] = 3, then knl[j] denotes the
       normal integrated sextupole strength scaled by Bρo.
"""
@inline function multipole_kick!(i, v, work, mm, knl, ksl)
  @assert size(work, 2) >= 3 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 3) for multipole_kick!()."
  @inbounds begin #@FastGTPSA! begin
  jm = length(mm)
  m = mm[jm]
  if m == 1 return v end
  work[i,1] = ar = knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
  work[i,2] = ai = knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
  jm -= 1
  while m > 2
    m -= 1
    work[i,3] = work[i,1] * v[i,XI] - work[i,2] * v[i,YI]
    work[i,2] = work[i,1] * v[i,YI] + work[i,2] * v[i,XI]
    work[i,1] = work[i,3]
    if m == mm[jm]
      work[i,1] += knl[jm] * v[i,XI] - ksl[jm] * v[i,YI]
      work[i,2] += knl[jm] * v[i,YI] + ksl[jm] * v[i,XI]
      jm -= 1
    end
  end
  v[i,PXI] -= work[i,1]
  v[i,PYI] += work[i,2]
  end #end
  return v
end # function multipole_kick!()


#function binom(m::Integer, x, y)
#  """
#  This function computes the real and imaginary parts of
#  (x + i y)^m. One can use these components to compute,
#  among other things, the multipole kick induced by two-
#  dimensional multipole magnets.
#  """
#  if m == 0
#    return [ 1.0, 0.0 ]
#  end
#  ar = x
#  ai = y
#  mm = m
#  while mm > 1
#    mm -= 1
#    t  = x * ar - y * ai
#    ai = y * ar + x * ai
#    ar = t
#  end
#  return [ ar, ai ]
#end # function binom()


#
# ===============  E X A C T   S E C T O R   B E N D  ===============
#
"""
This function implements exact symplectic tracking through a
sector bend, derived using the Hamiltonian (25.9) given in the
BMad manual. As a consequence of using that Hamiltonian, the
reference value of βγ must be that of a particle with the
design energy.  Should we wish to change that, we shall need
to carry both reference and design values.

Arguments
—————————
beta_0: β_0 = (βγ)_0 / √(γ_0^2)
brho_0: Bρ_0,  reference magnetic rigidity
hc: coordinate frame curvature
b0: magnet field strength
e1: entrance face angle (+ve angle <=> toward rbend)
e2: exit face angle (+ve angle <=> toward rbend)
Lr: element arc length
"""
@inline function exact_sbend!(i, v, work, beta_0, brho_0, hc, b0, e1, e2, Lr)
  @assert size(work, 2) >= 5 && size(work, 1) == size(v, 1) "Size of work matrix must be at least ($size(v, 1), 5) for multipole_kick!()."
  @inbounds begin @FastGTPSA! begin
  rho = brho0 / b0
  ang = hc * Lr
  c1 = cos(ang)
  s1 = sin(ang)

  work[i,1] = sqrt((1.0 + v[i,PZI])^2 - (v[i,PXI]^2 + v[i,PYI]^2))  # P_s
  work[i,2] = sqrt((1.0 + v[i,PZI])^2 - v[i,PYI]^2)                 # P_α
  work[i,3] = (1.0 + hc * v[i,XI]) / (hc * rho)                     # scaled (1 + h x)
  work[i,4] = work[i,1] - work[i,3]                                 # Px'/h
  work[i,5] = ang + asin(v[i,PXI] / work[i,2]) - asin((v[i,PXI] * c1 + work[i,4] * s1) / work[i,2])  # α + φ1 - φ2
  # high-precision computation of x-final
  v[i,XI] = (v[i,XI] * c1 - Lr * sin(ang / 2) * sincu(ang / 2)
             + rho * (v[i,PXI] + ((v[i,PXI]^2 + (work[i,1] + work[i,4]) * work[i,3]) * s1 - 2v[i,PXI] * work[i,4] * c1)
                             / (sqrt(work[i,2]^2 - (v[i,PXI] * c1 + work[i,4] * s1)^2) + work[i,1] * c1)) * s1)
  v[i,PXI] = v[i,PXI] * c1 + work[i,4] * s1
  v[i,YI] = v[i,YI] + rho * v.py * work[i,5]
  # high-precision computation of z-final
  v[i,ZI] = (v[i,ZI] - rho * (1.0 + v[i,PZI]) * work[i,5]
               + (1.0 + v[i,PZI]) * Lr / (beta_0 * sqrt(1.0 / beta_0^2 + (2 + v[i,PZI]) * v[i,PZI])))

  end end
  return v
end # function exact_sbend!()


end
