using ArgParse
using LinearAlgebra
using Plots
using Printf
using WriteVTK

include("Basic.jl")
include("ODESolver.jl")
include("Physical.jl")
include("WriteIO.jl")

#using .Basic
#using .ODESolver
#using .Physical

function parse_commandline()

  s = ArgParseSettings()
  @add_arg_table s begin
    "--nx", "-n"
    help = "number of points in each direction"
    arg_type = Int
    default = 101
    "--ngh"
    help = "number of ghost points in each direction"
    arg_type = Int
    default = 2
    "--out_every"
    help = "output every so many steps"
    arg_type = Int
    default = 10
    "--cfl"
    help = "Courant factor"
    arg_type = Float64
    default = 0.25
    "--itlast"
    help = "maximum time steps"
    arg_type = Int
    default = 200
    "--ggif"
    help = "if generate gifs"
    arg_type = Bool
    default = false
  end
  return parse_args(s)

end

function main()

  println("===================================================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("===================================================================")

  params = parse_commandline()
  nx = params["nx"]
  ngh = params["ngh"]
  itlast = params["itlast"]
  out_every = params["out_every"]
  ggif = params["ggif"]
  cfl = params["cfl"]

  bbox = [[-4.0, 4.0], [-1.0, 1.0]]
  nbuf = ngh * 4
  grid = Basic.Grid(nx, bbox, ngh, nbuf, cfl)
  gfs = Basic.GridFunction(2, grid)

  ###############
  # Intial Data #
  ###############
  println("Setting up initial conditions...")
  Physical.InitialData!(gfs)

  @printf("Simulation time: %.4f, iteration %d. E = %.4f\n",
          gfs.grid.time, 0, Physical.Energy(gfs))

  WriteIO.dump("data", gfs, 0)

  ##########
  # Evolve #
  ##########
  println("Start evolution...")
  for i = 1:itlast
    ODESolver.Evolve!(Physical.WaveRHS!, gfs)
    @printf("Simulation time: %.4f, iteration %d. E = %.4f\n",
            gfs.grid.time, i, Physical.Energy(gfs))

    if (mod(i, out_every) == 0)
      WriteIO.dump("data", gfs, i)
    end
  end

  ########
  # Exit #
  ########
  println("-------------------------------------------------------------------")
  println("  Successfully Done")
  println("-------------------------------------------------------------------")

  #################
  # Generate gifs #
  #################
  if ggif
    WriteIO.generate_gifs("data")
  end

end

main()
