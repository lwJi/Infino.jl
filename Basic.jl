module Basic

export Grid, GridF

struct Grid

  xmin :: Float64
  xmax :: Float64
  nx :: Int64
  x :: Vector{Float64}
  dx :: Float64

  function Grid(nx0, bbox)
    nx = nx0
    xmin = bbox[1]
    xmax = bbox[2]
    dx = (xmax - xmin) / (nx -1)
    x = LinRange(xmin, xmax, nx)
    new(xmin, xmax, nx, x, dx)
  end

end

struct GridF

  function GridF()
  end

end

end # module Basic
