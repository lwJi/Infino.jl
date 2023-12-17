module Sync

include("Algo.jl")

function Prolongation(gfs, l, interp_in_time::Bool; ord_s = 1, ord_t = 1)
    nxa = gfs.grid.levs[l].nxa
    nbuf = gfs.grid.levs[l].nbuf
    if2c = gfs.grid.levs[l].if2c
    aligned = gfs.grid.levs[l].aligned
    levsfs = gfs.levs

    for j = 1:2  # left or right
        for v = 1:gfs.nd
            uf = levsfs[l].u[v]
            uc_p = levsfs[l-1].u_p[v]
            if (interp_in_time)
                uc = levsfs[l-1].u[v]
                uc_pp = levsfs[l-1].u_pp[v]
                for i = 1:nbuf
                    f = (j == 1) ? i : nxa - i + 1
                    c = if2c[f]
                    ucs =
                        (aligned[f]) ? [uc_pp[c], uc_p[c], uc[c]] :
                        ucs = [
                            Algo.Interpolation(uc_pp, c, ord_s),
                            Algo.Interpolation(uc_p, c, ord_s),
                            Algo.Interpolation(uc, c, ord_s),
                        ]
                    uf[f] = Algo.Interpolation(ucs, 2, ord_t)
                end
            else
                # here we assume that we always march coarse level first: l in 1:lmax
                for i = 1:nbuf
                    f = (j == 1) ? i : nxa - i + 1
                    c = if2c[f]
                    uf[f] = ((aligned[f]) ? uc_p[c] : Algo.Interpolation(uc_p, c, ord_s))
                end
            end
        end
    end
end

function Restriction(gfs, l)
    nxa = gfs.grid.levs[l+1].nxa
    nbuf = gfs.grid.levs[l+1].nbuf
    if2c = gfs.grid.levs[l+1].if2c
    aligned = gfs.grid.levs[l+1].aligned
    levsfs = gfs.levs
    for v = 1:gfs.nd
        uf = levsfs[l+1].u[v]
        uc = levsfs[l].u[v]
        for f = 1+nbuf:nxa-nbuf
            if aligned[f]
                uc[if2c[f]] = uf[f]
            end
        end
    end
end

end
