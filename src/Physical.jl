module Physical

function InitialData!(gfs)

  println("Setting initial data ...")

  nx = gfs.grid.nx
  u = gfs.u[1]
  rho = gfs.u[2]
  x = gfs.grid.x

  amp = 1.0
  sig = 0.1
  x0 = 0.0

  for i = 1:nx
    u[i] = amp * exp(-((x[i] - x0) / sig)^2)
    rho[i] = 0.0
  end

end

function WaveScalarRHS!()
end

end # module Physical
