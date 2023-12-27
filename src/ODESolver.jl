module ODESolver

include("Sync.jl")

#===============================================================================
Eovlve!:
    * evolve one complete time step for all levels
===============================================================================#
function Evolve!(f::Function, gfs; Mongwane = false, apply_trans_zone = false)
    lmax = length(gfs.levs)

    #----------------------------------------#
    # march the first substep for all levels #
    #----------------------------------------#
    for l = 1:lmax  # notice that we march coarse level first
        if l > 1
            Mongwane ? Sync.Prolongation_Mongwane(gfs, l, false) :
            Sync.Prolongation(gfs, l, false)
            if apply_trans_zone
                Sync.ApplyTransitionZone(gfs, l, false)
            end
        end
        Mongwane ? rk4_Mongwane!(f, gfs.levs[l]) : rk4!(f, gfs.levs[l])
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
                        Sync.Restriction(gfs, l; apply_trans_zone)  # from l+1 to l
                    end
                    # from l-1 to l
                    Mongwane ?
                    Sync.Prolongation_Mongwane(gfs, l, mod(substeps[l], 2) == 0) :
                    Sync.Prolongation(gfs, l, mod(substeps[l], 2) == 0)
                    if apply_trans_zone
                        Sync.ApplyTransitionZone(gfs, l, mod(substeps[l], 2) == 0)
                    end
                    Mongwane ? rk4_Mongwane!(f, gfs.levs[l]) : rk4!(f, gfs.levs[l])
                end
            end
        end
    end

    #------------------------#
    # Restriction all levels #
    #------------------------#
    for l = lmax-1:-1:1  # notice that we restrict fine level first
        Sync.Restriction(gfs, l; apply_trans_zone)
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

    @. u_pp = u_p * 1.0
    @. u_p = u * 1.0
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

    @. u_pp = u_p * 1.0
    @. u_p = u * 1.0
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

function rk4_Mongwane!(f::Function, levf)
    u = levf.u
    u_p = levf.u_p
    r = levf.rhs
    w = levf.w
    k1 = levf.k[1]
    k2 = levf.k[2]
    k3 = levf.k[3]
    k4 = levf.k[4]
    lev = levf.lev
    t = lev.time
    dt = lev.dt
    isrt = lev.is_lev1 ? 1 : 1 + lev.nbuf
    iend = lev.is_lev1 ? lev.nxa : lev.nxa - lev.nbuf

    @. u_p = u * 1.0
    lev.time = t
    f(lev, r, u)
    for v = 1:levf.nd
        k1[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k1 * (1 / 6)

    @. w = u_p + k1 * (1 / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    for v = 1:levf.nd
        k2[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k2 * (1 / 3)

    @. w = u_p + k2 * (1 / 2)
    lev.time = t + 0.5 * dt
    f(lev, r, w)
    for v = 1:levf.nd
        k3[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k3 * (1 / 3)

    @. w = u_p + k3
    lev.time = t + dt
    f(lev, r, w)
    for v = 1:levf.nd
        k4[v][isrt:iend] = r[v][isrt:iend] * dt
    end
    @. u += k4 * (1 / 6)
    lev.time = t + dt
end

end
