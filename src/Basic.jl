module Basic

export Grid, GridF

struct Grid

    xmin::Float64
    xmax::Float64
    nx::Int64
    x::Vector{Float64}
    dx::Float64

    function Grid(nx0, bbox)
        nx = nx0
        xmin = bbox[1]
        xmax = bbox[2]
        dx = (xmax - xmin) / (nx - 1)
        x = LinRange(xmin, xmax, nx)
        new(xmin, xmax, nx, x, dx)
    end

end

struct GridF

    ndim::Int64
    grid::Grid
    u::Array{Array{Float64,1},1}
    dxu::Array{Array{Float64,1},1}

    function GridF(ndim, grid)

        nx = grid.nx
        u = Array{Array{Float64,1},1}(undef, ndim)
        dxu = Array{Array{Float64,1},1}(undef, ndim)
        for i = 1:ndim
            u[i] = zeros(Float64, nx)
            dxu[i] = zeros(Float64, nx)
        end
        new(ndim, grid, u, dxu)

    end

end

end # module Basic
