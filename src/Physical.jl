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

  @. psi = amp * exp(-((x - x0) / sig)^2)
  @. Pi = 0.0

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

  @. psi_rhs = Pi
  @. Pi_rhs = ddpsi

end

function Energy(gfs)

  nx = gfs.grid.nx
  dx = gfs.grid.dx
  psi = gfs.u[1]
  Pi = gfs.u[2]
  dpsi = zeros(Float64, nx)
  Derivs.derivs_1st!(dpsi, psi, dx, 4)

  E::Float64 = 0.0
  for i = 1:nx
    E += (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi[i] * dpsi[i])
  end

  return E * dx

end

end # module Physical
