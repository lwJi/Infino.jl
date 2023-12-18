module InitialData

include("Derivs.jl")
include("Sync.jl")
include("ODESolver.jl")
include("Physical.jl")

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
    for l = lmax-1:-1:1
        Sync.Restriction(gfs, l)
    end
end

function NegativeWaveRHS!(lev, r, u)
    Physical.WaveRHS!(lev, r, u)
    @. r = -r
end

function MarchBackwards!(gfs)
    for l = 1:length(gfs.levs)
        levf = gfs.levs[l]

        if l > 1
            Sync.Prolongation(gfs, l, false)
        end
        ODESolver.rk4!(NegativeWaveRHS!, levf)

        u = levf.u
        u_p = levf.u_p
        u_pp = levf.u_pp
        @. u_pp = u_p
        @. u_p = u
        @. u = u_pp

        levf.lev.time = 0.0
    end
    gfs.grid.time = 0.0
end

end # module InitialData
