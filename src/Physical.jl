module Physical

function InitialData!(gf)

    println("Setting initial data ...")

    nx = gf.grid.nx
    u = gf.u[1]
    rho = gf.u[2]
    x = gf.grid.x

    amp = 1.0
    sig = 1.0
    x0 = 0.0

    for i = 1:nx
        u[i] = amp * exp(-((x[i] - x0) / sig)^2)
        rho[i] = 0.0
    end
end

end # module Physical
