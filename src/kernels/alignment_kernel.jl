@makekernel fastgtpsa=true function track_alignment_straight_entering!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  track_translation!(i, coords, (x_off, y_off, 0.0), 0.0)
  q = inv_rot_quaternion(x_rot, y_rot, tilt)
  L2 = 0.5 * L * ele_orient
  dz_new = coord_rotation!(i, coords, q, -L2 - z_off)
  alignment_drift!(i, coords, -dz_new - L2)
end

#

@makekernel fastgtpsa=true function track_alignment_straight_exiting!(i, coords::Coords,
                                            x_off, y_off, z_off, x_rot, y_rot, tilt, ele_orient, L)
  v = coords.v
  L2 = 0.5 * L * ele_orient

  q = rot_quaternion(x_rot, y_rot, tilt)
  dz_new = coord_rotation!(i, coords, q, L2)

  alive = (coords.state[i] == STATE_ALIVE)
  v[i,XI] = vifelse(alive, v[i,XI] + x_off, v[i,XI])
  v[i,YI] = vifelse(alive, v[i,YI] + y_off, v[i,YI])

  alignment_drift!(i, coords, L2 - z_off - dz_new)
end

#---------------------------------------------------------------------------------------------------

"""
    track_alignment_bend_entering(i, coords::Coords, x_off, y_off, z_off, x_rot, y_rot, tilt, 
                                              g_ref, tilt_ref, ele_orient, L) -> dr, q

Returns `dr` origin shift and `q` quaternion rotation for the coordinate transformation
from the nominal bend entrance face (in branch coordinates) to the actual entrance face
(in body coordinates) taking into account the element alignment parameters.

## Arguments
- `x_off`, `y_off`, `z_off`   Element offset.
- `x_rot`, `y_rot`, 'tilt`    Element orientation.
- `g_ref`                     Reference g = 1/bend radius.
- `tilt_ref`                  Branch coords tilt.
- `ele_orient`                Element longitudinal orientation: +1 => normal, -1 => reversed.
- `L`                         Element length.

## Returns
- `dr`    Coordinate origin shift.
- `q`     Quaternion rotation.
"""

@makekernel fastgtpsa=true function track_alignment_bend_entering!(i, coords::Coords,
                          x_off, y_off, z_off, x_rot, y_rot, tilt, g_ref, tilt_ref, ele_orient, L)
  # Transform to coord system with origin at element chord midpoint.

  r, q = coord_bend_transform(i, coords, 0.5*L, g_ref, tilt_ref)
  println("***A: $r  : $q")

  dr = -L^2 * g_ref * one_cos_norm(tilt) * (cos(tilt_ref), sin(tilt_ref), 0.0)
  r += dr
  println("***B: $r  : $q")

  # Misalignment transform

  dr = (-x_off, -y_off, 0.0)
  r += dr
  println("***C: $r  : $q")
 
  dq = inv_rot_quaternion(x_rot, y_rot, tilt)
  q = quat_mul(q, dq)
  r = quat_rotate(r, dq)
  println("***C: $r  : $q")

  # Transform to body coord system with origin at entrance.

  dq = inv_rot_quaternion(0.0, 0.0, tilt_ref)
  q = quat_mul(q, dq)
  r = quat_rotate(r, dq)
  println("***D: $r  : $q")

  dr = L^2 * g_ref * one_cos_norm(tilt) * (1.0, 0.0, 0.0)
  r += dr
  println("***E: $r  : $q")


  r, q = coord_bend_transform(-0.5*L, g_ref, 0.0)
  println("***E: $r  : $q")

  # Transform coords by (r, q) and drift to entrance

  track_translation!(i, coords, r, 0.0)
  dz_new = coord_rotation!(i, coords, q, -r[3])
  println("***AA: $(v[1,:]) : $z1")

  alignment_drift!(i, coords, -dz_new)
  println("***BB: $(v[1,:]) : $z1")
end
 