module Sync

include("Algo.jl")
include("Symb.jl")

#===============================================================================
Prolongation_new: use Mongwane's method
    * from level l-1 to level l
    * we assume that we always march coarse level first (for l in 2:lmax)
===============================================================================#
function Prolongation_new(gfs, l, interp_in_time::Bool; ord_s = 3, ord_t = 2)
    levf = gfs.levs[l]
    levc = gfs.levs[l-1]

    for j = 1:2  # left or right
        for v = 1:gfs.nd
            for i = 1:nbuf
                f = (j == 1) ? i : nxa - i + 1
                c = if2c[f]
                if aligned[f]
                    kfs = calc_kfs_from_kcs(
                        [levc.k[m][v][c] for m = 1:4],
                        dtc,
                        interp_in_time,
                    )
                    for m = 1:3
                        levf.k[m][v][f] = kfs[m]
                    end
                else
                    kfss = zeros(Float64, 3, 4)
                    for k = 1:4
                        kfss[:, k] = calc_kfs_from_kcs(
                            [levc.k[m][v][c+k-2] for m = 1:4],
                            dtc,
                            interp_in_time,
                        )
                    end
                    for m = 1:3
                        levf.k[m][v][f] = Algo.Interpolation(kfss[m, :], 2, ord_s)
                    end
                end
            end
        end
    end
end

function calc_kfs_from_kcs(kcs, dtc, interp_in_time::Bool)
    t0_f = (interp_in_time) ? 0.5 : 0.0
    thalf_f = (interp_in_time) ? 0.75 : 0.25
    dtf = 0.5 * dtc
    d1yc_t0 = Symb.dy(1)(t0_f, dtc, kcs)
    d1yc = Symb.dy(1)(thalf_f, dtc, kcs)
    d2yc = Symb.dy(2)(thalf_f, dtc, kcs)
    d3yc = Symb.dy(3)(thalf_f, dtc, kcs)
    fyd2yc = 4 * (kcs[3] - kcs[2]) / dtc^3
    return [
        dtf * d1yc_t0,
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc - fyd2yc),
        dtf * d1yc + 0.5 * dtf^2 * d2yc + 0.125 * dtf^3 * (d3yc + fyd2yc),
    ]
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
                    ucs =
                        (aligned[f]) ? [uc_pp[c], uc_p[c], uc[c]] :
                        [
                            Algo.Interpolation(uc_pp, c, ord_s),
                            Algo.Interpolation(uc_p, c, ord_s),
                            Algo.Interpolation(uc, c, ord_s),
                        ]
                    uf[f] = Algo.Interpolation(ucs, 2, ord_t)
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
function Restriction(gfs, l)
    nxa = gfs.grid.levs[l+1].nxa
    nbuf = gfs.grid.levs[l+1].nbuf
    if2c = gfs.grid.levs[l+1].if2c
    aligned = gfs.grid.levs[l+1].aligned
    for v = 1:gfs.nd
        uf = gfs.levs[l+1].u[v]
        uc = gfs.levs[l].u[v]
        for f = 1+nbuf:nxa-nbuf  # only interior
            if aligned[f]
                uc[if2c[f]] = uf[f]
            end
        end
    end
end

end
