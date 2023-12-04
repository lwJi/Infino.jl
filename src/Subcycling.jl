using Printf
using WriteVTK
using Plots

include("Basic.jl")
include("ODESolver.jl")
include("Physical.jl")

using .Basic
using .ODESolver
using .Physical

function main()

  println("==================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("==================================")

  nx = 101
  bbox = [-1.0, 1.0]
  grid = Basic.Grid(nx, bbox)

  gfs = Basic.GridF(2, grid)
  println("nd = ", gfs.nd)
  println("nx = ", gfs.grid.nx)
  println("dx = ", gfs.grid.dx)
  println("dt = ", gfs.grid.dt)

  println("time = ", gfs.grid.time)

  ###############
  # Intial Data #
  ###############
  println("Setting Initial Data ...")
  Physical.InitialData!(gfs)

  # output
  a_u = Animation()
  plt_u = plot(gfs.grid.x, gfs.u[2], labal="u")
  # plt_u = scatter!(gfs.grid.x, gfs.u[1])
  frame(a_u, plt_u)

  ##########
  # Evolve #
  ##########
  println("Evolving ...")
  for i = 1:9
    ODESolver.euler!(Physical.WaveRHS!, gfs)

    println("time = ", gfs.grid.time)
    plt_u = plot(gfs.grid.x, gfs.u[2], labal="u")
    # plt_u = scatter!(gfs.grid.x, gfs.u[1])
    frame(a_u, plt_u)
  end

  # output
  gif(a_u, "u.gif")

  ########
  # Exit #
  ########
  println("Successfully Done ...")

end

main()
