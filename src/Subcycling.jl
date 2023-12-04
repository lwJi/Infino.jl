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
  bbox = [-1.0, 1.0]
  grid = Basic.Grid(nx, bbox)

  gfs = Basic.GridF(2, grid)
  println()
  println("  nd = ", gfs.nd)
  println("  nx = ", gfs.grid.nx)
  println("  dx = ", gfs.grid.dx)
  println("  dt = ", gfs.grid.dt)
  println()

  ###############
  # Intial Data #
  ###############
  Physical.InitialData!(gfs)

  @printf("Simulation time: %.4f, iteration %d. |psi| = %.4f\n",
          gfs.grid.time, 0, norm(gfs.u[1]))

  a_u = Animation()
  plt_u = plot(gfs.grid.x, gfs.u[1], labal="u")
  # plt_u = scatter!(gfs.grid.x, gfs.u[1])
  frame(a_u, plt_u)

  ##########
  # Evolve #
  ##########
  itlast = 20
  out_every = 20

  for i = 1:itlast
    ODESolver.rk4!(Physical.WaveRHS!, gfs)
    @printf("Simulation time: %.4f, iteration %d. |psi| = %.4f\n",
            gfs.grid.time, i, norm(gfs.u[1]))

    if (mod(i, out_every) == 0)
      plt_u = plot(gfs.grid.x, gfs.u[1], labal="u")
      # plt_u = scatter!(gfs.grid.x, gfs.u[1])
      frame(a_u, plt_u)
    end
  end

  # output
  gif(a_u, "u.gif")

  ########
  # Exit #
  ########
  println("-------------------------------------------------------------------")
  println("  Successfully Done")
  println("-------------------------------------------------------------------")

end

main()
