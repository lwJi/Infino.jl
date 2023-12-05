using Plots

include("../Basic.jl")

function main()

  nx = 101
  bnds = [[-4.0, 4.0], [0.5, 1.0]]
  grid = Basic.Grid(nx, bnds)

  println("lev1:")
  println("  nx   = ", grid.levs[1].nx)
  println("  xmin = ", grid.levs[1].xmin)
  println("  xmax = ", grid.levs[1].xmax)
  println("  dx   = ", grid.levs[1].dx)
  println("  dt   = ", grid.levs[1].dt)
  # println("  x    = ", grid.levs[1].x)
  println("lev2:")
  println("  nx   = ", grid.levs[2].nx)
  println("  xmin = ", grid.levs[2].xmin)
  println("  xmax = ", grid.levs[2].xmax)
  println("  dx   = ", grid.levs[2].dx)
  println("  dt   = ", grid.levs[2].dt)
  # println("  x    = ", grid.levs[2].x)

  ax = Animation()
  plt = plot(grid.levs[1].x, zeros(Float64, grid.levs[1].nx),
             xlim=(0.3, 1.2), ylim=(-0.1, 0.6))
  plt = scatter!(grid.levs[1].x, zeros(Float64, grid.levs[1].nx))
  plt = scatter!(grid.levs[2].x, zeros(Float64, grid.levs[2].nx) .+ 0.1)
  frame(ax, plt)

  gif(ax, "x.gif")

end

main()
