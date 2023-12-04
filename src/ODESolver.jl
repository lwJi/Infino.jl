module ODESolver

function Evolve!(gfs, itlast, out_every)
end

############################
# Time Integration Methods #
############################
function euler!(f::Function, gfs)

  u = gfs.u
  u_p = gfs.u_p
  r = gfs.rhs
  grid = gfs.grid
  t = grid.time
  dt = grid.dt
  dx = grid.dx
  x = grid.x

  u_p = u
  grid.time = t
  f(grid, r, u)
  @. u += r * dt
  grid.time = t + dt

end

function rk4!(f::Function, gfs)

  u = gfs.u
  u_p = gfs.u_p
  r = gfs.rhs
  w = gfs.w
  grid = gfs.grid
  t = grid.time
  dt = grid.dt
  dx = grid.dx
  x = grid.x

  u_p = u
  grid.time = t
  f(grid, r, u)
  @. u += r * (dt/6)

  @. w = u_p + r * (dt/2)
  grid.time = t + 0.5 * dt
  f(grid, r, w)
  @. u += r * (dt/3)

  @. w = u_p + r * (dt/2)
  grid.time = t + 0.5 * dt
  f(grid, r, w)
  @. u += r * (dt/3)

  @. w = u_p + r * dt
  grid.time = t + dt
  f(grid, r, w)
  @. u += r * (dt/6)
  grid.time = t + dt

end

end
