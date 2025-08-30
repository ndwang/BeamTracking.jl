@makekernel fastgtpsa=true function track_aperture_rectangular!(i, coords::Coords, x1, x2, y1, y2)
  v = coords.v

  coords.state[i] = vifelse(v[i,XI] < x1, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,XI] > x2, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,YI] < y1, STATE_LOST, coords.state[i])
  coords.state[i] = vifelse(v[i,YI] > y2, STATE_LOST, coords.state[i])
end

#

@makekernel fastgtpsa=true function track_aperture_elliptical!(i, coords::Coords, x1, x2, y1, y2)
  v = coords.v

  x0 = 0.5 * (x2 + x1) 
  y0 = 0.5 * (y2 + y1) 
  xw = 0.5 * (x2 - x1) 
  yw = 0.5 * (y2 - y1) 
  r = ((v[i,XI] - x0) / xw)^2 + ((v[i,YI] - y0) / yw)^2

  coords.state[i] = vifelse(r > 1 || xw < 0 || yw < 0, STATE_LOST, coords.state[i])
end
