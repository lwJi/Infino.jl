module ODESolver

function euler!(f::Function, gfs)
  nd = gfs.nd

  t = gfs.grid.time
  dt = gfs.grid.dt
  dx = gfs.grid.dx
  x = gfs.grid.x
  u = gfs.u
  u_p = gfs.u_p
  r = gfs.r

  u_p = u
  gfs.grid.time = t
  f(gfs, r, u)
  u = u + r * dt
  gfs.grid.time = t + dt

end

function rk4!()
end

end
