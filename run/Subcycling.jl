using Infino
using Printf
using TOML

function main(pars, out_dir)
    println(
        "--------------------------------------------------------------------------------",
    )
    println()
    println(" =======")
    println("         \\")
    println("                         *******************")
    println("          \\                Subcycling Test")
    println("  -----       =======    *******************")
    println("         \\ \\")
    println("           --  ----- ")
    println("         /")
    println(" =======   -- =======")
    println()
    println(
        "--------------------------------------------------------------------------------",
    )

    ########################
    # Read Parameter Files #
    ########################
    nx = pars["parameters"]["nx"]
    ngh = pars["parameters"]["ngh"]
    nbuf = pars["parameters"]["nbuf"]
    itlast = pars["parameters"]["itlast"]
    out_every = pars["parameters"]["out_every"]
    bbox = pars["parameters"]["bbox"]
    cfl = haskey(pars["parameters"], "cfl") ? pars["parameters"]["cfl"] : 0.25
    diss = haskey(pars["parameters"], "diss") ? pars["parameters"]["diss"] : 0.0
    subcycling =
        haskey(pars["parameters"], "subcycling") ? pars["parameters"]["subcycling"] : true
    Mongwane =
        haskey(pars["parameters"], "Mongwane") ? pars["parameters"]["Mongwane"] : false
    ntrans = haskey(pars["parameters"], "ntrans") ? pars["parameters"]["ntrans"] : 3
    apply_trans_zone =
        haskey(pars["parameters"], "apply_trans_zone") ?
        pars["parameters"]["apply_trans_zone"] : false
    initial_data =
        haskey(pars["parameters"], "initial_data") ? pars["parameters"]["initial_data"] :
        "Gaussian"
    println("Parameters:")
    println("  cfl        = ", cfl)
    println("  Mongwane   = ", Mongwane)
    println("  trans_zone = ", apply_trans_zone)
    println("  itlast     = ", itlast)
    println("  out_every  = ", out_every)
    println("  out_dir    = ", out_dir)

    ########################
    # build grid structure #
    ########################
    grid = Infino.Basic.Grid(
        nx,
        bbox,
        ngh,
        nbuf;
        ntrans = ntrans,
        cfl = cfl,
        diss = diss,
        subcycling = subcycling,
    )
    gfs = Infino.Basic.GridFunction(2, grid)

    ###############
    # Intial Data #
    ###############
    println("Setting up initial conditions...")
    println("  initial data type: $initial_data")
    if initial_data == "Gaussian"
        Infino.InitialData.Gaussian!(gfs)
    elseif initial_data == "sinusoidal"
        Infino.InitialData.sinusoidal!(gfs)
    else
        println("Initial data type '$initial_data' unsupported yet")
        exit()
    end

    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    if !Mongwane
        Infino.InitialData.MarchBackwards!(gfs)
        Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    end

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
        Infino.ODESolver.Evolve!(
            Infino.Physical.WaveRHS!,
            gfs;
            Mongwane = Mongwane,
            apply_trans_zone = apply_trans_zone,
        )

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
    # Done #
    ########
    println(
        "--------------------------------------------------------------------------------",
    )
    println("Successfully Done.")
end

function redirect_to_files(dofunc, outfile, errfile)
    open(outfile, "w") do out
        open(errfile, "w") do err
            redirect_stdout(out) do
                redirect_stderr(err) do
                    dofunc()
                end
            end
        end
    end
end

#===============================================================================
Start Execution
===============================================================================#
if length(ARGS) < 1
    println("Usage: julia Subcycling.jl parfile.toml")
    exit(1)
end
pars_path = ARGS[1]
pars = TOML.parsefile(pars_path)

# create output directory
out_dir = joinpath(
    dirname(pars_path),
    haskey(pars["parameters"], "out_dir") ? pars["parameters"]["out_dir"] :
    splitext(basename(pars_path))[1],
)
if isdir(out_dir)
    println("Removing old directory '$out_dir'...")
    rm(out_dir, recursive = true)
end
println("Creating new directory '$out_dir'...")
mkdir(out_dir)

# copy parfile into out_dir
cp(pars_path, out_dir * "/" * basename(pars_path))

# config
redirect_std =
    haskey(pars["configs"], "redirect_std") ? pars["configs"]["redirect_std"] : true

if redirect_std
    # redirect output and error
    redirect_to_files("./stdout.txt", "./stderr.txt") do
        main(pars, out_dir)
    end
    mv("./stdout.txt", out_dir * "/stdout.txt")
    mv("./stderr.txt", out_dir * "/stderr.txt")
else
    main(pars, out_dir)
end
#===============================================================================
End Execution
===============================================================================#
