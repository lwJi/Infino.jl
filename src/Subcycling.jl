include("Basic.jl")

using .Basic

function main()

    println("==================================")
    println("  Welcome to Subcycling Test !!!  ")
    println("==================================")

    nx = 21
    bbox = [-1.0, 1.0]
    grid = Basic.Grid(nx, bbox)

    println("nx = ", grid.nx)

end

main()
