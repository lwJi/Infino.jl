using Infino
using Printf
using TOML

function main()

    println("===================================================================")
    println("  Welcome to Subcycling Test !!!  ")
    println("===================================================================")

    if length(ARGS) < 1
        println("Usage: julia Subcycling.jl parfile.toml")
        exit(1)
    end

    ########################
    # Read Parameter Files #
    ########################
    pars_path = ARGS[1]
    pars = TOML.parsefile(pars_path)
    bbox = pars["parameters"]["bbox"]
    cfl = pars["parameters"]["cfl"]
    nx = pars["parameters"]["nx"]
    ngh = pars["parameters"]["ngh"]
    itlast = pars["parameters"]["itlast"]
    out_every = pars["parameters"]["out_every"]
    bbox = pars["parameters"]["bbox"]
    if haskey(pars["parameters"], "out_dir")
        out_dir_base = pars["parameters"]["out_dir"]
    else
        out_dir_base = splitext(basename(pars_path))[1]
    end
    out_dir = joinpath(dirname(pars_path), out_dir_base)
    # print pars
    println("Parameters:")
    println("  cfl       = ", cfl)
    println("  itlast    = ", itlast)
    println("  out_every = ", out_every)
    println("  out_dir   = ", out_dir)
    # create output dir
    if isdir(out_dir)
        println("Removing old directory '$out_dir'...")
        rm(out_dir, recursive = true)
    end
    println("Creating new directory '$out_dir'...")
    mkdir(out_dir)
    # build grid structure
    nbuf = ngh * 4
    grid = Infino.Basic.Grid(nx, bbox, ngh, nbuf; cfl = cfl)
    gfs = Infino.Basic.GridFunction(2, grid)

    ###############
    # Intial Data #
    ###############
    println("Setting up initial conditions...")
    Infino.InitialData.Gaussian!(gfs)
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    Infino.InitialData.MarchBackwards!(gfs)
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    @printf(
        "Simulation time: %.4f, iteration %d. E = %.4f\n",
        gfs.grid.time,
        0,
        Infino.Physical.Energy(gfs)
    )

    Infino.WriteIO.dump(out_dir, gfs, 0)

    ##########
    # Evolve #
    ##########
    println("Start evolution...")
    for i = 1:itlast
        Infino.ODESolver.Evolve!(Infino.Physical.WaveRHS!, gfs)
        Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
        @printf(
            "Simulation time: %.4f, iteration %d. E = %.4f\n",
            gfs.grid.time,
            i,
            Infino.Physical.Energy(gfs)
        )

        if (mod(i, out_every) == 0)
            Infino.WriteIO.dump(out_dir, gfs, i)
        end
    end

    ########
    # Exit #
    ########
    println("-------------------------------------------------------------------")
    println("  Successfully Done")
    println("-------------------------------------------------------------------")

end

main()
