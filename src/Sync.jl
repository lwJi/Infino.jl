module Sync

include("Algo.jl")
include("DenseOutput.jl")

#===============================================================================
Functions needed by Mongwane's subcycling method
===============================================================================#
function calc_kfs_from_kcs(kcs, dtc, interp_in_time::Bool)
    t0_f = interp_in_time ? 0.5 : 0.0
    dtf = 0.5 * dtc
    d1yc = DenseOutput.dy1(t0_f, dtc, kcs)
    d2yc = DenseOutput.dy2(t0_f, dtc, kcs)
    d3yc = DenseOutput.dy3(t0_f, dtc, kcs)
    fyd2yc = 4 * (kcs[3] - kcs[2]) / dtc^3
    return [
        dtf * d1yc,
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc - fyd2yc),
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc + fyd2yc),
    ]
end

function transition_profile(a, b, x; type = 1)
    @assert(b - x > -1e-12, "b less then a is not allowed")
    @assert(x - a > -1e-12, "b less then a is not allowed")
    t0 = (x - a) / (b - a)
    t = t0 < 0 ? 0 : (t0 > 1 ? 1 : t0)

    if type == 1  # boxstep
        return t
    elseif type == 2  # smoothstep
        return 3 * t^2 - 2 * t^3
    elseif type == 3  # smootherstep
        return 10 * t^3 - 15 * t^4 + 6 * t^5
    else
        println("Transition profile (type $type) is not supported yet.")
        exit()
    end
end

#===============================================================================
ApplyTransitionZone: apply transition zone
===============================================================================#
function ApplyTransitionZone(gfs, l, interp_in_time::Bool; ord_s = 3)
    nxa = gfs.grid.levs[l].nxa
    nbuf = gfs.grid.levs[l].nbuf
    ntrans = gfs.grid.levs[l].ntrans
    if2c = gfs.grid.levs[l].if2c
    aligned = gfs.grid.levs[l].aligned
    levf = gfs.levs[l]
    levc = gfs.levs[l-1]
    # for transition zone
    xbox = gfs.grid.levs[l].xbox
    dxf = gfs.grid.levs[l].dx

    for j = 1:2  # left or right
        a = (j == 1) ? xbox[1] : xbox[2] - (ntrans - 1) * dxf
        b = (j == 1) ? xbox[1] + (ntrans - 1) * dxf : xbox[2]
        for v = 1:gfs.nd
            uf = levf.u[v]
            uc_p = levc.u_p[v]
            for i = 1:ntrans
                f = (j == 1) ? i + nbuf : nxa - i + 1 - nbuf
                c = if2c[f]
                w = transition_profile(a, b, gfs.levs[l].x[f])
                if aligned[f]
                    kcs = [levc.k[m][v][c] for m = 1:4]
                    ys = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                    uf[f] = (1 - w) * ys + w * uf[f]
                else
                    ys = zeros(Float64, 4)
                    for ic = 1:4
                        kcs = [levc.k[m][v][c+ic-2] for m = 1:4]
                        ys[ic] =
                            interp_in_time ? DenseOutput.y(0.5, uc_p[c+ic-2], kcs) :
                            uc_p[c+ic-2]
                    end
                    uf[f] = (1 - w) * Algo.Interpolation(ys, 2, ord_s) + w * uf[f]
                end
            end
        end
    end
end

#===============================================================================
Prolongation_Mongwane: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function Prolongation_Mongwane(gfs, l, interp_in_time::Bool; ord_s = 3)
    nxa = gfs.grid.levs[l].nxa
    nbuf = gfs.grid.levs[l].nbuf
    if2c = gfs.grid.levs[l].if2c
    aligned = gfs.grid.levs[l].aligned
    dtc = gfs.grid.levs[l-1].dt
    levf = gfs.levs[l]
    levc = gfs.levs[l-1]

    for j = 1:2  # left or right
        for v = 1:gfs.nd
            uf = levf.u[v]
            uc_p = levc.u_p[v]
            for i = 1:nbuf
                f = (j == 1) ? i : nxa - i + 1
                c = if2c[f]
                if aligned[f]
                    kcs = [levc.k[m][v][c] for m = 1:4]
                    kfs = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                    # setting k
                    for m = 1:3
                        levf.k[m][v][f] = kfs[m]
                    end
                    # setting u
                    uf[f] = interp_in_time ? DenseOutput.y(0.5, uc_p[c], kcs) : uc_p[c]
                else
                    kfss = zeros(Float64, 3, 4)
                    ys = zeros(Float64, 4)
                    for ic = 1:4
                        kcs = [levc.k[m][v][c+ic-2] for m = 1:4]
                        kfss[:, ic] = calc_kfs_from_kcs(kcs, dtc, interp_in_time)
                        ys[ic] =
                            interp_in_time ? DenseOutput.y(0.5, uc_p[c+ic-2], kcs) :
                            uc_p[c+ic-2]
                    end
                    # setting k
                    for m = 1:3
                        levf.k[m][v][f] = Algo.Interpolation(kfss[m, :], 2, ord_s)
                    end
                    # setting u
                    uf[f] = Algo.Interpolation(ys, 2, ord_s)
                end
            end
        end
    end
end

#===============================================================================
Prolongation:
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function Prolongation(gfs, l, interp_in_time::Bool; ord_s = 3, ord_t = 2)
    nxa = gfs.grid.levs[l].nxa
    nbuf = gfs.grid.levs[l].nbuf
    if2c = gfs.grid.levs[l].if2c
    aligned = gfs.grid.levs[l].aligned
    levf = gfs.levs[l]
    levc = gfs.levs[l-1]

    for j = 1:2  # left or right
        for v = 1:gfs.nd
            uf = levf.u[v]
            uc_p = levc.u_p[v]
            if interp_in_time
                uc = levc.u[v]
                uc_pp = levc.u_pp[v]
                for i = 1:nbuf
                    f = (j == 1) ? i : nxa - i + 1
                    c = if2c[f]
                    if aligned[f]
                        uf[f] = Algo.Interpolation([uc_pp[c], uc_p[c], uc[c]], 2, ord_t)
                    else
                        ucss = zeros(Float64, 3, 4)
                        for ic = 1:4
                            ucss[:, ic] = [uc_pp[c+ic-2], uc_p[c+ic-2], uc[c+ic-2]]
                        end
                        uf[f] = Algo.Interpolation(
                            [Algo.Interpolation(ucss[m, :], 2, ord_s) for m = 1:3],
                            2,
                            ord_t,
                        )
                    end
                end
            else
                for i = 1:nbuf
                    f = (j == 1) ? i : nxa - i + 1
                    c = if2c[f]
                    uf[f] = ((aligned[f]) ? uc_p[c] : Algo.Interpolation(uc_p, c, ord_s))
                end
            end
        end
    end
end

#===============================================================================
Restriction:
    * from level l+1 to level l
    * we assume that we always march fine level first (for l in lmax-1:-1:1)
    * we assume all the levels are at the same time slice
===============================================================================#
function Restriction(gfs, l; apply_trans_zone = false)
    nxa = gfs.grid.levs[l+1].nxa
    nbuf = gfs.grid.levs[l+1].nbuf
    ntrans = gfs.grid.levs[l+1].ntrans
    if2c = gfs.grid.levs[l+1].if2c
    aligned = gfs.grid.levs[l+1].aligned
    isrt = apply_trans_zone ? 1 + nbuf + ntrans : 1 + nbuf
    iend = apply_trans_zone ? nxa - nbuf - ntrans : nxa - nbuf
    for v = 1:gfs.nd
        uf = gfs.levs[l+1].u[v]
        uc = gfs.levs[l].u[v]
        for f = isrt:iend  # only interior
            if aligned[f]
                uc[if2c[f]] = uf[f]
            end
        end
    end
end

end
