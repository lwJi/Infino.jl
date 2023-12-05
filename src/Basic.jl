module Basic

export Level, LevelFunction, Grid, GridFunction

struct Level

  nx::Int64
  xmin::Float64
  xmax::Float64
  dx::Float64
  dt::Float64

  function Level(nx, bnd, dt)
    xmin = bnd[1]
    xmax = bnd[2]
    dx = (xmax - xmin) / (nx - 1)
    new(nx, xmin, xmax, dx, dt)
  end

end

struct LevelFunction
  nd ::Int64
  lev::Level
  x  ::Array{Float64,1}
  u  ::Array{Array{Float64,1},1}
  u_p::Array{Array{Float64,1},1}
  rhs::Array{Array{Float64,1},1}
  w  ::Array{Array{Float64,1},1}

  function LevelFunction(nd, lev)
    x   = LinRange(lev.xmin, lev.xmax, lev.nx)
    u   = Array{Array{Float64,1},1}(undef, nd)
    u_p = Array{Array{Float64,1},1}(undef, nd)
    rhs = Array{Array{Float64,1},1}(undef, nd)
    w   = Array{Array{Float64,1},1}(undef, nd)

    nx = lev.nx
    for i = 1:nd
      u[i]   = zeros(Float64, nx)
      u_p[i] = zeros(Float64, nx)
      rhs[i] = zeros(Float64, nx)
      w[i]   = zeros(Float64, nx)
    end
    new(nd, lev, x, u, u_p, rhs, w)
  end
end

mutable struct Grid

  levs::Vector{Level}
  time::Float64
  dt  ::Float64

  function Grid(nx1, bnds::Vector{Vector{Float64}}, cfl=0.4, t=0.0)
    dx1 = (bnds[1][2] - bnds[1][1]) / (nx1 - 1)
    dt1 = cfl * dx1
    lev1 = Level(nx1, bnds[1], dt1)

    levs = Vector{Level}([lev1])
    for i = 2:length(bnds)
      dx = dx1 / 2^(i-1)
      levl = levs[i-1]
      xnew = LinRange(levl.xmin, levl.xmax, (levl.nx - 1)*2 + 1)
      xmin = xnew[argmin(abs.(xnew .- bnds[i][1]))]
      ncell = floor(Int, (bnds[i][2] - xmin) / dx)
      bnd = [xmin, xmin + ncell * dx]
      push!(levs, Level(ncell+1, bnd, cfl * dx))
    end

    println("Grid Structure:")
    for i = 1:length(levs)
      println("lev[", i, "],")
      println("  nx   = ", levs[i].nx)
      println("  xmin = ", levs[i].xmin)
      println("  xmax = ", levs[i].xmax)
      println("  dx   = ", levs[i].dx)
      println("  dt   = ", levs[i].dt)
    end

    new(levs, t, dt1)
  end

end

struct GridFunction

  nd   ::Int64
  grid ::Grid
  levfs::Vector{LevelFunction}

  function GridFunction(nd, grid)
    levfs = [LevelFunction(nd, lev) for lev in grid.lev]
    new(nd, grid, levfs)
  end

end

end # module Basic
