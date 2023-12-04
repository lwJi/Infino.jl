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
    "--out_every"
      help = "output every so many steps"
      arg_type = Int
      default = 20
    "--cfl"
      help = "Courant factor"
      arg_type = Float64
      default = 0.25
    "--itlast"
      help = "maximum time steps"
      arg_type = Int
      default = 500
  end
  return parse_args(s)

end

function main()

  println("===================================================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("===================================================================")

  params = parse_commandline()
  nx = params["nx"]
  itlast = params["itlast"]
  out_every = params["out_every"]
  cfl = params["cfl"]

  bbox = [-4.0, 4.0]
  grid = Basic.Grid(nx, bbox, cfl)

  gfs = Basic.GridF(2, grid)
  println()
  println("Grid structure and so on:")
  println("  nd = ", gfs.nd)
  println("  nx = ", gfs.grid.nx)
  println("  dx = ", gfs.grid.dx)
  println("  dt = ", gfs.grid.dt)
  println()

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

  plt_psi = plot(gfs.grid.x, gfs.u[1], ylim=(-1,1), labal="psi")
  # plt_psi = scatter!(gfs.grid.x, gfs.u[1])
  plt_Pi = plot(gfs.grid.x, gfs.u[2], ylim=(-4,4), labal="Pi")
  frame(a_psi, plt_psi)
  frame(a_Pi, plt_Pi)

  ##########
  # Evolve #
  ##########
  println("Start evolution...")
  for i = 1:itlast
    ODESolver.rk4!(Physical.WaveRHS!, gfs)
    @printf("Simulation time: %.4f, iteration %d. E = %.4f\n",
            gfs.grid.time, i, Physical.Energy(gfs))

    if (mod(i, out_every) == 0)
      plt_psi = plot(gfs.grid.x, gfs.u[1], ylim=(-1,1), labal="psi")
      # plt_psi = scatter!(gfs.grid.x, gfs.u[1])
      plt_Pi = plot(gfs.grid.x, gfs.u[2], ylim=(-4,4), labal="Pi")
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
