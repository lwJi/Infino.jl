module Basic

export Level, LevelFunction, Grid, GridFunction

mutable struct Level

  nx  ::Int64
  ngh ::Int64
  nxa ::Int64
  ibox::Array{Int64,1}
  xbox::Array{Float64,1}
  dx  ::Float64
  dt  ::Float64
  time::Float64

  function Level(nx, ngh, xbox, dt, t)
    nxa = nx + 2*ngh
    ibox = [1 + ngh, nxa - ngh]
    dx = (xbox[2] - xbox[1]) / (nx - 1)
    new(nx, ngh, nxa, ibox, xbox, dx, dt, t)
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
    xmin = lev.xbox[1] - lev.ngh * lev.dx
    xmax = lev.xbox[2] + lev.ngh * lev.dx
    x   = LinRange(xmin, xmax, lev.nxa)
    u   = Array{Array{Float64,1},1}(undef, nd)
    u_p = Array{Array{Float64,1},1}(undef, nd)
    rhs = Array{Array{Float64,1},1}(undef, nd)
    w   = Array{Array{Float64,1},1}(undef, nd)
    for i = 1:nd
      u[i]   = zeros(Float64, lev.nxa)
      u_p[i] = zeros(Float64, lev.nxa)
      rhs[i] = zeros(Float64, lev.nxa)
      w[i]   = zeros(Float64, lev.nxa)
    end
    new(nd, lev, x, u, u_p, rhs, w)
  end

end

mutable struct Grid

  levs::Vector{Level}
  dt  ::Float64
  time::Float64

  function Grid(nx1, ngh, xboxs::Vector{Vector{Float64}}, cfl=0.4, t=0.0)
    dx1 = (xboxs[1][2] - xboxs[1][1]) / (nx1 - 1)
    dt1 = cfl * dx1
    lev1 = Level(nx1, ngh, xboxs[1], dt1, t)
    levs = Vector{Level}([lev1])
    for i = 2:length(xboxs)
      dx = dx1 / 2^(i-1)
      levl = levs[i-1]
      xnew = LinRange(levl.xbox[1], levl.xbox[2], (levl.nx - 1)*2 + 1)
      xmin = xnew[argmin(abs.(xnew .- xboxs[i][1]))]
      ncell = floor(Int, (xboxs[i][2] - xmin) / dx)
      xbox = [xmin, xmin + ncell * dx]
      push!(levs, Level(ncell+1, ngh, xbox, cfl * dx, t))
    end
    println("Grid Structure:")
    for i = 1:length(levs)
      println("lev[", i, "],")
      println("  nx   = ", levs[i].nx)
      println("  xmin = ", levs[i].xbox[1])
      println("  xmax = ", levs[i].xbox[2])
      println("  dx   = ", levs[i].dx)
      println("  dt   = ", levs[i].dt)
    end
    new(levs, dt1, t)
  end

end

struct GridFunction

  nd  ::Int64
  grid::Grid
  levs::Vector{LevelFunction}

  function GridFunction(nd, grid)
    levfs = [LevelFunction(nd, lev) for lev in grid.levs]
    new(nd, grid, levfs)
  end

end

end # module Basic
