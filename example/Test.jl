using Plots

include("../Basic.jl")

function main()

    nx = 101
    ngh = 2
    nbuf = ngh * 3
    bnds = [[-4.0, 4.0], [0.6, 1.4], [0.8, 1.2]]
    grid = Basic.Grid(nx, bnds, ngh, nbuf)

    ax = Animation()
    x1 = LinRange(grid.levs[1].xbox[1], grid.levs[1].xbox[2], grid.levs[1].nx)
    x2 = LinRange(grid.levs[2].xbox[1], grid.levs[2].xbox[2], grid.levs[2].nx)
    x3 = LinRange(grid.levs[3].xbox[1], grid.levs[3].xbox[2], grid.levs[3].nx)
    plt = plot(x1, zeros(Float64, grid.levs[1].nx), xlim = (0.3, 1.7), ylim = (-0.1, 0.6))
    plt = scatter!(x1, zeros(Float64, grid.levs[1].nx))
    plt = scatter!(x2, zeros(Float64, grid.levs[2].nx) .+ 0.1)
    plt = scatter!(x3, zeros(Float64, grid.levs[3].nx) .+ 0.2)
    frame(ax, plt)

    gif(ax, "x.gif")

end

main()
