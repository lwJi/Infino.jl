module Physical

include("Derivs.jl")
using .Derivs

function InitialData!(gfs)

  nx = gfs.grid.nx
  psi = gfs.u[1]
  Pi = gfs.u[2]
  x = gfs.grid.x

  amp = 1.0
  sig = 0.1
  x0 = 0.0

  for i = 1:nx
    psi[i] = amp * exp(-((x[i] - x0) / sig)^2)
    Pi[i] = 0.0
  end

end

function WaveRHS!(grid, r, u)

  nx = grid.nx
  dx = grid.dx
  psi = u[1]
  Pi = u[2]
  psi_rhs = r[1]
  Pi_rhs = r[2]

  # derivatives
  ddpsi = zeros(Float64, nx)
  Derivs.derivs_2nd!(ddpsi, psi, dx, 4)

  # dtu   = rho
  # dtrho = ddu
  for i = 1:nx
    psi_rhs[i] = Pi[i]
    Pi_rhs[i] = ddpsi[i]
  end

end

end # module Physical
