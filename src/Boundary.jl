module Boundary

function ApplyPeriodicBoundaryCondition!(gfs)
    lev1fs = gfs.levs[1]
    nxa = gfs.grid.levs[1].nxa
    nbuf = gfs.grid.levs[1].nbuf
    for v = 1:gfs.nd
        u = lev1fs.u[v]
        for i = 1:nbuf
            u[i] = u[nxa-2*nbuf+i]
        end
        for i = nxa:-1:nxa-nbuf+1
            u[i] = u[2*nbuf-nxa+i]
        end
    end
end

end
