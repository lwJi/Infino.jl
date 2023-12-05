using Plots

include("../Basic.jl")

function main()

  nx = 101
  bnds = [[-4.0, 4.0], [0.5, 1.0]]
  grid = Basic.Grid(nx, bnds)

  ax = Animation()
  plt = plot(grid.levs[1].x, zeros(Float64, grid.levs[1].nx),
             xlim=(0.3, 1.2), ylim=(-0.1, 0.6))
  plt = scatter!(grid.levs[1].x, zeros(Float64, grid.levs[1].nx))
  plt = scatter!(grid.levs[2].x, zeros(Float64, grid.levs[2].nx) .+ 0.1)
  frame(ax, plt)

  gif(ax, "x.gif")

end

main()
