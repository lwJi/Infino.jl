module Basic

export Level, LevelFunction, Grid, GridFunction

mutable struct Level

    nx::Int64  # num of interior grid points
    ngh::Int64  # num of ghost points
    nbuf::Int64  # num of buffer points
    ntrans::Int64  # num of transition zone points
    nxa::Int64  # num of all grid points
    fdord::Int64  # finite difference order
    xbox::Array{Float64,1}  # size computational domain (interior)
    dx::Float64
    dt::Float64
    time::Float64
    diss::Float64
    is_lev1::Bool
    # dimension nxa
    if2c::Array{Int64,1}  # map between indexes of current and its parent level
    aligned::Array{Bool,1}  # if grid aligned with coarse grid

    function Level(nx, ngh, nbuf, ntrans, fdord, xbox, dt, t, diss, is_lev1, if2c, aligned)
        nxa = nx + 2 * nbuf
        dx = (xbox[2] - xbox[1]) / (nx - 1)
        new(
            nx,
            ngh,
            nbuf,
            ntrans,
            nxa,
            fdord,
            xbox,
            dx,
            dt,
            t,
            diss,
            is_lev1,
            if2c,
            aligned,
        )
    end

end

struct LevelFunction

    nd::Int64  # num of evolution variables
    lev::Level  # level structure
    x::Array{Float64,1}  # coordinates
    u::Array{Array{Float64,1},1}  # state vectors
    u_p::Array{Array{Float64,1},1}  # previous state vectors
    u_pp::Array{Array{Float64,1},1}  # previous previous state vectors
    rhs::Array{Array{Float64,1},1}  # rhs of state vectors
    w::Array{Array{Float64,1},1}  # intermediate state vectors
    # intermediate state vectors for new subcycling
    k::Array{Array{Array{Float64,1},1},1}

    function LevelFunction(nd, lev)
        noffset = (lev.nxa - lev.nx) / 2  # take account of buffer zone
        xmin = lev.xbox[1] - noffset * lev.dx
        xmax = lev.xbox[2] + noffset * lev.dx
        x = LinRange(xmin, xmax, lev.nxa)
        u = Array{Array{Float64,1},1}(undef, nd)
        u_p = Array{Array{Float64,1},1}(undef, nd)
        u_pp = Array{Array{Float64,1},1}(undef, nd)
        rhs = Array{Array{Float64,1},1}(undef, nd)
        w = Array{Array{Float64,1},1}(undef, nd)
        k = Array{Array{Array{Float64,1},1},1}(undef, 4)
        for j = 1:4
            k[j] = Array{Array{Float64,1},1}(undef, nd)
        end

        for i = 1:nd
            u[i] = zeros(Float64, lev.nxa)
            u_p[i] = zeros(Float64, lev.nxa)
            u_pp[i] = zeros(Float64, lev.nxa)
            rhs[i] = zeros(Float64, lev.nxa)
            w[i] = zeros(Float64, lev.nxa)
            for j = 1:4
                k[j][i] = zeros(Float64, lev.nxa)
            end
        end
        new(nd, lev, x, u, u_p, u_pp, rhs, w, k)
    end

end

mutable struct Grid

    levs::Vector{Level}
    dt::Float64
    time::Float64
    subcycling::Bool  # turn on subcycling or not

    function Grid(
        nx1,  # num of interior grid points at base level
        xboxs::Vector{Vector{Float64}},
        ngh,
        nbuf;
        ntrans = 0,
        fdord = 4,
        cfl = 0.4,
        t = 0.0,
        diss = 0.0,
        subcycling = true,
        verbose = true,
    )
        # build the first level (base level)
        dx1 = (xboxs[1][2] - xboxs[1][1]) / (nx1 - 1)
        dt1 = subcycling ? cfl * dx1 : cfl * dx1 / 2^(length(xboxs) - 1)
        lev1 = Level(nx1, ngh, nbuf, ntrans, fdord, xboxs[1], dt1, t, diss, true, [], [])
        levs = Vector{Level}([lev1])
        # build the rest levels
        for i = 2:length(xboxs)
            dx = dx1 / 2^(i - 1)
            dt = (subcycling ? cfl * dx : dt1)
            levl = levs[i-1]  # level lower than the current level (parent level)
            # if we refine parent level everywhere
            xl = LinRange(levl.xbox[1], levl.xbox[2], (levl.nx - 1) * 2 + 1)
            # find those two which are closest to current level boundaries
            imin = argmin(abs.(xl .- xboxs[i][1]))
            imax = argmin(abs.(xl .- xboxs[i][2]))
            #imin = findall(x->abs(x - xboxs[i][1]) <= dx + 1e-12, xl)[1]
            #imax = findall(x->abs(x - xboxs[i][2]) <= dx + 1e-12, xl)[end]
            xbox = [xl[imin], xl[imax]]
            nx = (imax - imin) + 1  # (floor(Int, (xl[imax] - xl[imin]) / dx)) + 1
            # maps between two levels
            if2c = div.(((imin-nbuf:imax+nbuf) .+ 1), 2) .+ nbuf
            aligned = mod.(((imin-nbuf:imax+nbuf) .+ 1), 2) .== 0
            # build level
            push!(
                levs,
                Level(
                    nx,
                    ngh,
                    nbuf,
                    ntrans,
                    fdord,
                    xbox,
                    dt,
                    t,
                    diss,
                    false,
                    if2c,
                    aligned,
                ),
            )
        end
        if verbose  # print
            println("Grid Structure:")
            for i = 1:length(levs)
                println("lev[", i, "],")
                println("  nx    = ", levs[i].nx)
                println("  ngh   = ", levs[i].ngh)
                println("  nbuf  = ", levs[i].nbuf)
                println("  fdord = ", levs[i].fdord)
                if length(levs[i].if2c) == levs[i].nxa
                    println(
                        "  ibox  = ",
                        [
                            levs[i].if2c[1+levs[i].nbuf],
                            levs[i].if2c[levs[i].nxa-levs[i].nbuf],
                        ],
                        ", ",
                        [
                            levs[i].aligned[1+levs[i].nbuf],
                            levs[i].aligned[levs[i].nxa-levs[i].nbuf],
                        ],
                    )
                end
                println("  xbox  = ", levs[i].xbox)
                println("  dx    = ", levs[i].dx)
                println("  dt    = ", levs[i].dt)
            end
        end
        # construct
        new(levs, dt1, t, subcycling)
    end

end

struct GridFunction

    nd::Int64
    grid::Grid
    levs::Vector{LevelFunction}

    function GridFunction(nd, grid)
        levfs = [LevelFunction(nd, lev) for lev in grid.levs]
        new(nd, grid, levfs)
    end

end

end # module Basic
