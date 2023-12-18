module InitialData

include("Derivs.jl")
include("Sync.jl")
include("ODESolver.jl")

function Gaussian!(gfs)
    lmax = length(gfs.levs)
    amp = 1.0
    sig = 0.2
    x0 = 0.0

    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x

        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end

    # for consistence
    for l = 2:lmax
        Sync.Prolongation(gfs, l, false)
    end
    for l = lmax-1:-1:1
        Sync.Restriction(gfs, l)
    end
end

function NegativeWaveRHS!(lev, r, u)
    psi = u[1]
    Pi = u[2]
    psi_rhs = r[1]
    Pi_rhs = r[2]
    # derivatives
    ddpsi = zeros(Float64, lev.nxa)
    Derivs.derivs_2nd!(ddpsi, psi, lev.dx, lev.ngh * 2)

    @. psi_rhs = -Pi
    @. Pi_rhs = -ddpsi
end

function MarchBackwards!(gfs)
    for l = 1:length(gfs.levs)
        levfs = gfs.levs[l]

        if l > 1
            Sync.Prolongation(gfs, l, false)
        end
        ODESolver.rk4!(NegativeWaveRHS!, levfs)

        u = levfs.u
        u_p = levfs.u_p
        u_pp = levfs.u_pp
        @. u_pp = u_p
        @. u_p = u
        @. u = u_pp
        levfs.lev.time = 0.0
    end
    gfs.grid.time = 0.0
end

end # module InitialData
