module Basic

export Grid, GridF

struct Grid

    nx::Int64
    xmin::Float64
    xmax::Float64
    dx::Float64
    x::Vector{Float64}

    function Grid(nx0, bbox)
        nx = nx0
        xmin = bbox[1]
        xmax = bbox[2]
        dx = (xmax - xmin) / (nx - 1)
        x = LinRange(xmin, xmax, nx)
        new(nx, xmin, xmax, dx, x)
    end

end

struct GridF

    nd::Int64
    grid::Grid
    u::Array{Array{Float64,1},1}
    dxu::Array{Array{Float64,1},1}

    function GridF(nd, grid)

        nx = grid.nx
        u = Array{Array{Float64,1},1}(undef, nd)
        dxu = Array{Array{Float64,1},1}(undef, nd)
        for i = 1:nd
            u[i] = zeros(Float64, nx)
            dxu[i] = zeros(Float64, nx)
        end
        new(nd, grid, u, dxu)

    end

end

end # module Basic
