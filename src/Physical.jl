module Physical

include("Derivs.jl")

#===============================================================================
WaveRHS!:
    * rhs of wave equation
        dot(psi) = Pi
        dot(Pi)  = ddpsi
===============================================================================#
function WaveRHS!(lev, r, u; interior_only = false)
    psi = u[1]
    Pi = u[2]
    psi_rhs = r[1]
    Pi_rhs = r[2]

    ddpsi = zeros(Float64, lev.nxa)
    psi_diss = zeros(Float64, lev.nxa)
    Pi_diss = zeros(Float64, lev.nxa)
    Derivs.derivs_2nd!(ddpsi, psi, lev.dx, lev.fdord)
    Derivs.derivs_diss!(psi_diss, psi, lev.dx, lev.fdord)
    Derivs.derivs_diss!(Pi_diss, Pi, lev.dx, lev.fdord)

    isrt = interior_only ? 1 + lev.nbuf : 1
    iend = interior_only ? lev.nxa - lev.nbuf : lev.nxa
    psi_rhs[isrt:iend] = Pi[isrt:iend] + lev.diss * psi_diss[isrt:iend]
    Pi_rhs[isrt:iend] = ddpsi[isrt:iend] + lev.diss * Pi_diss[isrt:iend]
end

#===============================================================================
Energy:
    * int_xmin^xmax (Pi^2/2 + dpsi^2/2)
    * calculate on base level (interior) only
===============================================================================#
function Energy(gfs)
    nxa = gfs.grid.levs[1].nxa
    nbuf = gfs.grid.levs[1].nbuf
    dx = gfs.grid.levs[1].dx
    psi = gfs.levs[1].u[1]
    Pi = gfs.levs[1].u[2]

    dpsi = zeros(Float64, nxa)
    Derivs.derivs_1st!(dpsi, psi, dx, gfs.grid.levs[1].fdord)

    E::Float64 = 0.0
    for i = 1+nbuf:nxa-nbuf
        E += (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi[i] * dpsi[i])
    end
    return E * dx
end

end # module Physical
