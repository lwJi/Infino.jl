using ArgParse
using LinearAlgebra
using Plots
using Printf
using WriteVTK

include("Basic.jl")
include("ODESolver.jl")
include("Physical.jl")

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
    "--xrange"
    help = "xlim for plots"
    arg_type = Tuple{Float64, Float64}
    default = (-4.7, 4.7)
    "--cfl"
    help = "Courant factor"
    arg_type = Float64
    default = 0.25
    "--itlast"
    help = "maximum time steps"
    arg_type = Int
    default = 200
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
  xrange = params["xrange"]
  cfl = params["cfl"]

  bbox = [[-4.0, 4.0], [-1.0, 1.0]]
  nbuf = ngh * 3
  grid = Basic.Grid(nx, bbox, ngh, nbuf, cfl)
  gfs = Basic.GridFunction(2, grid)

  yrange = (0, 1)
  a_psi = Animation()
  a_Pi = Animation()

  ###############
  # Intial Data #
  ###############
  println("Setting up initial conditions...")
  Physical.InitialData!(gfs)

  @printf("Simulation time: %.4f, iteration %d. E = %.4f\n",
          gfs.grid.time, 0, Physical.Energy(gfs))

  plt_psi = plot(gfs.levs[1].x, gfs.levs[1].u[1],
                 xlim=xrange, ylim=(-1,1), label="psi")
  for l = 1:length(gfs.levs)
    plt_psi = scatter!(gfs.levs[l].x, gfs.levs[l].u[1], label="")
  end
  plt_Pi = plot(gfs.levs[1].x, gfs.levs[1].u[2],
                xlim=xrange, ylim=(-4,4), label="Pi")
  frame(a_psi, plt_psi)
  frame(a_Pi, plt_Pi)

  ##########
  # Evolve #
  ##########
  println("Start evolution...")
  for i = 1:itlast
    ODESolver.Evolve!(Physical.WaveRHS!, gfs)
    @printf("Simulation time: %.4f, iteration %d. E = %.4f\n",
            gfs.grid.time, i, Physical.Energy(gfs))

    if (mod(i, out_every) == 0)
      plt_psi = plot(gfs.levs[1].x, gfs.levs[1].u[1],
                     xlim=xrange, ylim=(-1,1), label="psi")
      for l = 1:length(gfs.levs)
        plt_psi = scatter!(gfs.levs[l].x, gfs.levs[l].u[1], label="")
      end
      plt_Pi = plot(gfs.levs[1].x, gfs.levs[1].u[2],
                    xlim=xrange, ylim=(-4,4), label="Pi")
      frame(a_psi, plt_psi)
      frame(a_Pi, plt_Pi)
    end
  end

  # output
  gif(a_psi, "psi.gif")
  gif(a_Pi, "Pi.gif")

  ########
  # Exit #
  ########
  println("-------------------------------------------------------------------")
  println("  Successfully Done")
  println("-------------------------------------------------------------------")

end

main()
