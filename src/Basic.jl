module Basic

export Grid, GridF

struct Grid

  nx::Int64
  xmin::Float64
  xmax::Float64
  dx::Float64
  x::Vector{Float64}
  time::Float64
  dt::Float64

  function Grid(nx0, bbox, cfl=0.4, t=0.0)
    nx = nx0
    xmin = bbox[1]
    xmax = bbox[2]
    dx = (xmax - xmin) / (nx - 1)
    x = LinRange(xmin, xmax, nx)
    dt = dx * cfl
    new(nx, xmin, xmax, dx, x, t, dt)
  end

end

struct GridF

  nd::Int64
  grid::Grid
  u::Array{Array{Float64,1},1}
  u_p::Array{Array{Float64,1},1}
  r::Array{Array{Float64,1},1}
  w::Array{Array{Float64,1},1}

  function GridF(nd, grid)
    nx = grid.nx
    u = Array{Array{Float64,1},1}(undef, nd)
    u_p = Array{Array{Float64,1},1}(undef, nd)
    r = Array{Array{Float64,1},1}(undef, nd)
    w = Array{Array{Float64,1},1}(undef, nd)
    for i = 1:nd
      u[i] = zeros(Float64, nx)
      u_p[i] = zeros(Float64, nx)
      r[i] = zeros(Float64, nx)
      w[i] = zeros(Float64, nx)
    end
    new(nd, grid, u, u_p, r, w)
  end

end

end # module Basic
