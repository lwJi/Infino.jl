module ODESolver

include("Sync.jl")

function Evolve!(f::Function, gfs)
    tiny = 1e-12
    lmax = length(gfs.levs)

    #----------------------------------------#
    # march the first substep for all levels #
    #----------------------------------------#
    for l = 1:lmax
        if l > 1
            Sync.Prolongation(gfs, l, false)
        end
        rk4!(f, gfs.levs[l])
    end

    #-------------------------#
    # march the rest substeps #
    #-------------------------#
    substeps = ones(Int64, lmax)
    dt_min = gfs.grid.levs[lmax].dt
    for s = 2:2^(lmax-1)
        # march levels except coarest and finest ones
        for l = 2:lmax-1
            if (
                (abs(gfs.grid.levs[l+1].time - gfs.grid.levs[l].time) < tiny) &&
                (abs(gfs.grid.levs[1].time - gfs.grid.levs[l].time) > dt_min)
            )
                substeps[l] += 1
                Sync.Restriction(gfs, l)
                Sync.Prolongation(gfs, l, mod(substeps[l], 2) == 0)
                rk4!(f, gfs.levs[l])
            end
        end
        # march the finest level
        substeps[lmax] += 1
        Sync.Prolongation(gfs, lmax, mod(s, 2) == 0)
        rk4!(f, gfs.levs[lmax])
    end

    #------------------------#
    # Restriction all levels #
    #------------------------#
    for l = lmax-1:-1:1
        Sync.Restriction(gfs, l)
    end

    #------------------#
    # update grid time #
    #------------------#
    gfs.grid.time = gfs.grid.levs[1].time
end

############################
# Time Integration Methods #
############################
function euler!(f::Function, levfs)
    lev = levfs.lev
    u = levfs.u
    u_p = levfs.u_p
    r = levfs.rhs
    t = lev.time
    dt = lev.dt

    @. u_p = u
    lev.time = t
    f(lev, r, u)
    @. u += r * dt
    lev.time = t + dt
end

function rk4!(f::Function, levfs)
    lev = levfs.lev
    u = levfs.u
    u_p = levfs.u_p
    r = levfs.rhs
    w = levfs.w
    t = lev.time
    dt = lev.dt

    @. u_p = u
    lev.time = t
    f(lev, r, u)
    @. u += r * (dt / 6)

    @. w = u_p + r * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    @. u += r * (dt / 3)

    @. w = u_p + r * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    @. u += r * (dt / 3)

    @. w = u_p + r * dt
    lev.time = t + dt
    f(lev, r, w)
    @. u += r * (dt / 6)
    lev.time = t + dt
end

end
