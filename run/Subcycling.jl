using Infino
using Printf
using TOML

function main(pars, out_dir)

    println("===================================================================")
    println("  Welcome to Subcycling Test !!!  ")
    println("===================================================================")

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
    initial_data =
        haskey(pars["parameters"], "initial_data") ? pars["parameters"]["initial_data"] :
        "Gaussian"
    # print pars
    println("Parameters:")
    println("  cfl       = ", cfl)
    println("  itlast    = ", itlast)
    println("  out_every = ", out_every)
    println("  out_dir   = ", out_dir)
    # build grid structure
    nbuf = ngh * 4
    grid = Infino.Basic.Grid(nx, bbox, ngh, nbuf; cfl = cfl)
    gfs = Infino.Basic.GridFunction(2, grid)

    ###############
    # Intial Data #
    ###############
    println("Setting up initial conditions...")
    if initial_data == "Gaussian"
        Infino.InitialData.Gaussian!(gfs)
    elseif initial_data == "sinusoidal"
        Infino.InitialData.sinusoidal!(gfs)
    else
        println("Initial data type '$initial_data' unsupported yet")
        exit()
    end
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

###################
# Start Execution #
###################

if length(ARGS) < 1
    println("Usage: julia Subcycling.jl parfile.toml")
    exit(1)
end
pars_path = ARGS[1]
pars = TOML.parsefile(pars_path)

# create output directory
out_dir_base =
    haskey(pars["parameters"], "out_dir") ? pars["parameters"]["out_dir"] :
    splitext(basename(pars_path))[1]
out_dir = joinpath(dirname(pars_path), out_dir_base)
if isdir(out_dir)
    println("Removing old directory '$out_dir'...")
    rm(out_dir, recursive = true)
end
println("Creating new directory '$out_dir'...")
mkdir(out_dir)

# copy parfile
cp(pars_path, out_dir * "/" * basename(pars_path))

# redirect output and error
redirect_to_files("./stdout.txt", "./stderr.txt") do
    main(pars, out_dir)
end

mv("./stdout.txt", out_dir * "/stdout.txt")
mv("./stderr.txt", out_dir * "/stderr.txt")

#######
# End #
#######
