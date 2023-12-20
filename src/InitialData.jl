module InitialData

include("Derivs.jl")
include("Sync.jl")
include("ODESolver.jl")
include("Physical.jl")

#===============================================================================
Initial Data Types:
    * Gaussian
===============================================================================#
function Gaussian!(gfs; amp = 1.0, sig = 0.2, x0 = 0.0)
    lmax = length(gfs.levs)
    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x
        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        Sync.Restriction(gfs, l)
    end
end

function sinusoidal!(gfs)
    lmax = length(gfs.levs)
    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x
        @. psi = sin(2 * pi * (x - 0.0))
        @. Pi = -2 * pi * cos(2 * pi * (x - 0.0))
    end
    # restriction for consistence
    for l = lmax-1:-1:1
        Sync.Restriction(gfs, l)
    end
end

#===============================================================================
Spectial Treatment for Prolongation
    * evolve backwards to file u_p
===============================================================================#
function NegativeWaveRHS!(lev, r, u)
    Physical.WaveRHS!(lev, r, u)
    @. r = -r
end

function MarchBackwards!(gfs)
    for l = 1:length(gfs.levs)
        if l > 1
            Sync.Prolongation(gfs, l, false)
        end
        ODESolver.rk4!(NegativeWaveRHS!, gfs.levs[l])
        # save new u(-dt) -> u_p, u(0) -> u
        u = gfs.levs[l].u
        u_p = gfs.levs[l].u_p
        u_pp = gfs.levs[l].u_pp
        @. u_pp = u_p
        @. u_p = u
        @. u = u_pp
        gfs.levs[l].lev.time = 0.0
    end
    gfs.grid.time = 0.0
end

end # module InitialData
