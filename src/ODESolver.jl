module ODESolver

function Evolve!(f::Function, gfs)
  for l in 1:length(gfs.levs)
    rk4!(f, gfs.levs[l])
  end
end

############################
# Time Integration Methods #
############################
function euler!(f::Function, levfs)
  lev = levfs.lev
  u   = levfs.u
  u_p = levfs.u_p
  r   = levfs.rhs
  t  = lev.time
  dt = lev.dt

  u_p = u
  lev.time = t
  f(lev, r, u)
  @. u += r * dt
  lev.time = t + dt
end

function rk4!(f::Function, levfs)
  lev = levfs.lev
  u   = levfs.u
  u_p = levfs.u_p
  r   = levfs.rhs
  w   = levfs.w
  t  = lev.time
  dt = lev.dt

  u_p = u
  lev.time = t
  f(lev, r, u)
  @. u += r * (dt/6)

  @. w = u_p + r * (dt/2)
  lev.time = t + 0.5 * dt
  f(lev, r, w)
  @. u += r * (dt/3)

  @. w = u_p + r * (dt/2)
  lev.time = t + 0.5 * dt
  f(lev, r, w)
  @. u += r * (dt/3)

  @. w = u_p + r * dt
  lev.time = t + dt
  f(lev, r, w)
  @. u += r * (dt/6)
  lev.time = t + dt
end

end
