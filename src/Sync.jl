module Sync

function Prolongation(gfs, l, interp_in_time::Bool)
  nxa     = gfs.grid.levs[l].nxa
  nbuf    = gfs.grid.levs[l].nbuf
  if2c    = gfs.grid.levs[l].if2c
  aligned = gfs.grid.levs[l].aligned
  levsfs  = gfs.levs

  for v in 1:gfs.nd
    uf   = levsfs[l].u[v]
    uc_p = levsfs[l-1].u_p[v]
    if (interp_in_time)
      uc = levsfs[l-1].u[v]
      for f in 1:nbuf
        c = if2c[f]
        uf[f] = ((aligned[f])
                 ? (uc[c] + uc_p[c]) * 0.5
                 : (uc[c] + uc[c+1] + uc_p[c] + uc_p[c+1]) * 0.25)
      end
      for f in nxa:-1:nxa-nbuf+1
        c = if2c[f]
        uf[f] = ((aligned[f])
                 ? (uc[c] + uc_p[c]) * 0.5
                 : (uc[c] + uc[c+1] + uc_p[c] + uc_p[c+1]) * 0.25)
      end
    else
      # here we assume that we always march coarse level first: l in 1:lmax
      for f in 1:nbuf
        c = if2c[f]
        uf[f] = ((aligned[f]) ? uc_p[c] : (uc_p[c] + uc_p[c+1]) * 0.5)
      end
      for f in nxa:-1:nxa-nbuf+1
        c = if2c[f]
        uf[f] = ((aligned[f]) ? uc_p[c] : (uc_p[c] + uc_p[c+1]) * 0.5)
      end
    end
  end
end

function Restriction()
end

end
