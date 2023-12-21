module Sync

include("Algo.jl")

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
