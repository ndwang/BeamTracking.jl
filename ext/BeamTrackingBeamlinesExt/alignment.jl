# entering = 1 -> entering element.
# entering = -1 -> exiting element.

@inline function alignment(tm, bunch, alignmentP, bendP, L, entering::Int)
#  if !isactive(alignmentP); return nothing; end
#
#  x_off = alignmentP.x_off
#  y_off = alignmentP.y_off
#  z_off = alignmentP.z_off
#  x_rot = alignmentP.x_rot
#  y_rot = alignmentP.y_rot
#  tilt  = alignmentP.tilt
#  tilde_m, _, beta_0 = ExactTracking.drift_params(bunch.species, bunch.R_ref)
#  ele_orient = 1   # Need to extend this to reversed elements.
#
#  if isactive(bendP) && bendP.g_ref != 0
#    return KernelCall(track_alignment_bend!, (entering, x_off, y_off, z_off, 
#                x_rot, y_rot, tilt, bendP.g_ref, bendP.tilt_ref, ele_orient, tilde_m, beta_0, L))
#  else
#    return KernelCall(track_alignment_straight!, (entering, x_off, y_off, z_off, 
#                                         x_rot, y_rot, tilt, ele_orient, tilde_m, beta_0, L))
#  end
end