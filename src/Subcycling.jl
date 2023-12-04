using Printf
using WriteVTK
using Plots

include("Basic.jl")
include("Physical.jl")

using .Basic
using .Physical

function main()

  println("==================================")
  println("  Welcome to Subcycling Test !!!  ")
  println("==================================")

  nx = 101
  bbox = [-1.0, 1.0]
  grid = Basic.Grid(nx, bbox)
  println("nx = ", grid.nx)

  gfs = Basic.GridF(2, grid)
  println("nd = ", gfs.nd)

  Physical.InitialData!(gfs)

  time = [0.0]
  nfiles::Int64 = 0
  au = Animation()

  plt = plot(gfs.grid.x, gfs.u[1], labal="u")
  plt = scatter!(gfs.grid.x, gfs.u[1])
  frame(au, plt)

  gif(au, "u.gif")

  #fname = @sprintf("scalar_%04d", nfiles)
  #vtk_grid(fname, gfs.grid.x, LinRange(0.0, 0.0, 1)) do vtk
  #    vtk["u", VTKPointData()] = gfs.u[1]
  #    vtk["rho", VTKPointData()] = gfs.u[2]
  #    vtk["t", VTKPointData()] = time[1]
  #end

end

main()
