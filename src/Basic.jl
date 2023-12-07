module Basic

export Level, LevelFunction, Grid, GridFunction

mutable struct Level

  nx  ::Int64  # num of interior grid points
  ngh ::Int64  # num of ghost points
  nbuf::Int64  # num of buffer points
  nxa ::Int64  # num of all grid points
  xbox::Array{Float64,1}  # size computational domain (interior)
  dx  ::Float64
  dt  ::Float64
  time::Float64
  # dimension nxa
  ir2c::Array{Int64,1}  # map between indexes of current and its parent level
  align::Array{Bool,1}  # if grid align with coarse grid

  function Level(nx, ngh, nbuf, xbox, dt, t, ir2c, align)
    nxa = nx + 2*nbuf
    dx = (xbox[2] - xbox[1]) / (nx - 1)
    new(nx, ngh, nbuf, nxa, xbox, dx, dt, t, ir2c, align)
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
    noffset = (lev.nxa - lev.nx) / 2
    xmin = lev.xbox[1] - noffset * lev.dx
    xmax = lev.xbox[2] + noffset * lev.dx
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

  function Grid(nx1, xboxs::Vector{Vector{Float64}}, ngh, nbuf, cfl=0.4, t=0.0)
    dx1 = (xboxs[1][2] - xboxs[1][1]) / (nx1 - 1)
    dt1 = cfl * dx1
    lev1 = Level(nx1, ngh, nbuf, xboxs[1], dt1, t, [], [])
    levs = Vector{Level}([lev1])
    for i = 2:length(xboxs)
      dx = dx1 / 2^(i-1)
      levl = levs[i-1]
      xl = LinRange(levl.xbox[1], levl.xbox[2], (levl.nx - 1)*2 + 1)
      imin = argmin(abs.(xl .- xboxs[i][1]))
      imax = argmin(abs.(xl .- xboxs[i][2]))
      # ncell = floor(Int, (xboxs[i][2] - xmin) / dx)
      # xbox = [xmin, xmin + ncell * dx]
      # xl = LinRange(levl.xbox[1], levl.xbox[2], levl.nx)
      # imin = findall(x->abs(x - xboxs[i][1]) <= dx + 1e-12, xl)[1]
      # imax = findall(x->abs(x - xboxs[i][2]) <= dx + 1e-12, xl)[end]
      xbox = [xl[imin], xl[imax]]
      nx = (imax - imin) + 1  #  (floor(Int, (xl[imax] - xl[imin]) / dx)) + 1
      ir2c = div.(((imin-nbuf:imax+nbuf) .+ 1), 2)
      align = mod.(((imin-nbuf:imax+nbuf) .+ 1), 2) .== 0
      push!(levs, Level(nx, ngh, nbuf, xbox, cfl * dx, t, ir2c, align))
    end
    println("Grid Structure:")
    for i = 1:length(levs)
      println("lev[", i, "],")
      println("  nx   = ", levs[i].nx)
      println("  ngh  = ", levs[i].ngh)
      println("  nbuf = ", levs[i].nbuf)
      if (length(levs[i].ir2c) == levs[i].nxa)
        println("  ibox = ",
                [levs[i].ir2c[1+levs[i].nbuf],
                 levs[i].ir2c[levs[i].nxa-levs[i].nbuf]],
                ", ",
                [levs[i].align[1+levs[i].nbuf],
                 levs[i].align[levs[i].nxa-levs[i].nbuf]])
      end
      println("  xbox = ", levs[i].xbox)
      println("  dx   = ", levs[i].dx)
      println("  dt   = ", levs[i].dt)
      # println("  ir2c = ", levs[i].ir2c)
      # println("  align= ", levs[i].align)
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
