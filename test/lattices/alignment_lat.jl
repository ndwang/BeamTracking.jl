using BeamTracking
using Beamlines

@eles begin
  drift1 = Drift(L = 2, x_offset = 0.1, y_offset = 0.2, z_offset = 0.3, 
                    x_rot = 0.04, y_rot = 0.05, tilt = 0.06, tracking_method = Exact())

  drift2 = Drift(L = 2, tracking_method = Exact())

  bend1 = SBend(L = 2, g = 0.1, x_offset = 0.1, y_offset = 0.2, z_offset = 0.3, 
                    x_rot = 0.04, y_rot = 0.05, tilt = 0.06, tilt_ref = 0.5, tracking_method = Exact())
  bend2 = SBend(L = 2.0, g = 0.1, tilt_ref = 0., tracking_method = Exact())
end

bline_drift1 = Beamline([drift1], R_ref = 1.0, species_ref=Species("electron"))
bline_drift2 = Beamline([drift2], R_ref = 1.0, species_ref=Species("electron"))
bline_bend1 = Beamline([bend1], R_ref = 1.0, species_ref = Species("electron"))
bline_bend2 = Beamline([bend2], R_ref = 1.0, species_ref = Species("electron"))
;
