using Printf
using WriteVTK
using Plots
using LinearAlgebra

include("Basic.jl")
include("ODESolver.jl")
include("Physical.jl")

#using .Basic
#using .ODESolver
#using .Physical

function main()

  println("===================================================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("===================================================================")

  nx = 101
  bbox = [-4.0, 4.0]
  grid = Basic.Grid(nx, bbox)

  gfs = Basic.GridF(2, grid)
  println()
  println("Grid structure and so on:")
  println("  nd = ", gfs.nd)
  println("  nx = ", gfs.grid.nx)
  println("  dx = ", gfs.grid.dx)
  println("  dt = ", gfs.grid.dt)
  println()

  itlast = 500
  out_every = 20
  yrange = (0, 1)
  a_psi = Animation()
  a_Pi = Animation()

  ###############
  # Intial Data #
  ###############
  println("Setting up initial conditions...")
  Physical.InitialData!(gfs)

  @printf("Simulation time: %.4f, iteration %d. |psi| = %.4f\n",
          gfs.grid.time, 0, norm(gfs.u[1]))

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
    @printf("Simulation time: %.4f, iteration %d. |psi| = %.4f\n",
            gfs.grid.time, i, norm(gfs.u[1]))

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
