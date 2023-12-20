module ODESolver

include("Sync.jl")

#===============================================================================
Eovlve!:
    * evolve one complete time step for all levels
===============================================================================#
function Evolve!(f::Function, gfs)
    lmax = length(gfs.levs)

    #----------------------------------------#
    # march the first substep for all levels #
    #----------------------------------------#
    for l = 1:lmax  # notice that we march coarse level first
        if l > 1
            Sync.Prolongation(gfs, l, false)
        end
        rk4!(f, gfs.levs[l])
    end

    #-------------------------------------------------#
    # march the other substeps to the same time slice #
    #-------------------------------------------------#
    if gfs.grid.subcycling
        levs = gfs.grid.levs
        tiny = 1e-12
        dt_min = levs[lmax].dt
        substeps = ones(Int64, lmax)
        for s = 2:2^(lmax-1)  # from second to final substep (of the finest level)
            for l = 2:lmax  # march all levels except the coarest (from coarse to fine)
                if l == lmax || (
                    abs(levs[l].time - levs[l+1].time) < tiny &&
                    abs(levs[l].time - levs[1].time) > dt_min
                )
                    substeps[l] += 1
                    if l < lmax
                        Sync.Restriction(gfs, l)  # from l+1 to l
                    end
                    Sync.Prolongation(gfs, l, mod(substeps[l], 2) == 0)  # from l-1 to l
                    rk4!(f, gfs.levs[l])
                end
            end
        end
    end

    #------------------------#
    # Restriction all levels #
    #------------------------#
    for l = lmax-1:-1:1  # notice that we restrict fine level first
        Sync.Restriction(gfs, l)
    end

    #------------------#
    # update grid time #
    #------------------#
    gfs.grid.time = gfs.grid.levs[1].time
end

#===============================================================================
Time Integration Methods:
    * euler!, rk4!
===============================================================================#
function euler!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    u_pp = levf.u_pp
    r = levf.rhs
    lev = levf.lev
    t = lev.time
    dt = lev.dt

    @. u_pp = u_p
    @. u_p = u
    lev.time = t
    f(lev, r, u)
    @. u += r * dt
    lev.time = t + dt
end

function rk4!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    u_pp = levf.u_pp
    r = levf.rhs
    w = levf.w
    lev = levf.lev
    t = lev.time
    dt = lev.dt

    @. u_pp = u_p
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

function rk4_new!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    w = levf.w
    k1 = levf.k[1]
    k2 = levf.k[2]
    k3 = levf.k[3]
    k4 = levf.k[4]
    lev = levf.lev
    t = lev.time
    dt = lev.dt

    @. u_p = u
    lev.time = t
    f(lev, k1, u)
    @. u += k1 * (dt / 6)

    @. w = u_p + k1 * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, k2, w)
    @. u += k2 * (dt / 3)

    @. w = u_p + k2 * (dt / 2)
    lev.time = t + 0.5 * dt
    f(lev, k3, w)
    @. u += k3 * (dt / 3)

    @. w = u_p + k3 * dt
    lev.time = t + dt
    f(lev, k4, w)
    @. u += k4 * (dt / 6)
    lev.time = t + dt
end

end
