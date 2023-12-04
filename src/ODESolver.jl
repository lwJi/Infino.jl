module ODESolver

function Evolve!(gfs, itlast, out_every)
end

############################
# Time Integration Methods #
############################
function euler!(f::Function, gfs)

  u = gfs.u
  u_p = gfs.u_p
  r = gfs.r
  grid = gfs.grid
  t = grid.time
  dt = grid.dt
  dx = grid.dx
  x = grid.x

  u_p = u
  gfs.grid.time = t
  f(gfs.grid, r, u)
  @. u += r * dt
  gfs.grid.time = t + dt

end

function rk4!()
end

end
