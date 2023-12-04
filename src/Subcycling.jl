using Printf
using WriteVTK
using Plots

include("Basic.jl")
include("Derivs.jl")
include("Physical.jl")

using .Basic
using .Derivs
using .Physical

function main()

  println("==================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("==================================")

  nx = 101
  bbox = [-1.0, 1.0]
  grid = Basic.Grid(nx, bbox)

  gfs = Basic.GridF(2, grid)
  println("nd = ", gfs.nd)
  println("nx = ", gfs.grid.nx)
  println("dx = ", gfs.grid.dx)
  println("dt = ", gfs.grid.dt)

  println("time = ", gfs.grid.time)

  # intial data
  Physical.InitialData!(gfs)

  # derivatives
  du = Array{Array{Float64,1},1}(undef, gfs.nd)
  ddu = Array{Array{Float64,1},1}(undef, gfs.nd)
  for i = 1:gfs.nd
    du[i] = zeros(Float64, gfs.grid.nx)
    ddu[i] = zeros(Float64, gfs.grid.nx)
  end
  Derivs.calc_du!(du, gfs.u, gfs.grid.dx, 4)

  # output
  time = [0.0]
  nfiles::Int64 = 0
  a_u = Animation()
  plt_u = plot(gfs.grid.x, gfs.u[1], labal="u")
  plt_u = scatter!(gfs.grid.x, gfs.u[1])
  frame(a_u, plt_u)
  gif(a_u, "u.gif")
  a_du = Animation()
  plt_du = plot(gfs.grid.x, du[1], labal="du")
  plt_du = scatter!(gfs.grid.x, du[1])
  frame(a_du, plt_du)
  gif(a_du, "du.gif")

end

main()
