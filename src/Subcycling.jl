include("Basic.jl")
include("Physical.jl")

using .Basic
using .Physical

function main()

    println("==================================")
    println("  Welcome to Subcycling Test !!!  ")
    println("==================================")

    nx = 21
    bbox = [-1.0, 1.0]
    grid = Basic.Grid(nx, bbox)
    println("nx = ", grid.nx)

    gridf = Basic.GridF(2, grid)
    println("ndim = ", gridf.ndim)

    Physical.InitialData!(gridf)
end

main()
