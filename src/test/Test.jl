using Plots

include("../Basic.jl")

function main()

  nx = 101
  bnds = [[-4.0, 4.0], [0.6, 1.4], [0.8, 1.2]]
  grid = Basic.Grid(nx, bnds)

  ax = Animation()
  plt = plot(grid.xs[1], zeros(Float64, grid.levs[1].nx),
             xlim=(0.3, 1.7), ylim=(-0.1, 0.6))
  plt = scatter!(grid.xs[1], zeros(Float64, grid.levs[1].nx))
  plt = scatter!(grid.xs[2], zeros(Float64, grid.levs[2].nx) .+ 0.1)
  plt = scatter!(grid.xs[3], zeros(Float64, grid.levs[3].nx) .+ 0.2)
  frame(ax, plt)

  gif(ax, "x.gif")

end

main()
