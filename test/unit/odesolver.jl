using Infino

function example_ode!(lev, r, u)
    t = lev.time
    y = u[1]
    y_rhs = r[1]
    @. y_rhs = -0.5 * y + 0.5 * (98 + t) / (100 + t)^2
end

function analytical_solution(t)
    return (t + 100)^(-1) + 2 * exp(-0.5 * t)
end

@testset "Euler Method Tests" begin
    g = Infino.Basic.Grid(2, [[-1.0, 1.0]], 0, 0; cfl = 0.005, verbose = false)
    gfs = Infino.Basic.GridFunction(1, g)
    # initial data
    u = gfs.levs[1].u[1]
    x = gfs.levs[1].x
    @. u = analytical_solution.(0)
    # evolution
    for i = 1:4
        Infino.ODESolver.euler!(example_ode!, gfs.levs[1])
    end
    t = g.levs[1].time
    @test isapprox(gfs.levs[1].u[1][1], analytical_solution(t); rtol = 1e-4)
end

@testset "RK4 Method Tests" begin
    g = Infino.Basic.Grid(2, [[-1.0, 1.0]], 0, 0; cfl = 0.25, verbose = false)
    gfs = Infino.Basic.GridFunction(1, g)
    # initial data
    u = gfs.levs[1].u[1]
    x = gfs.levs[1].x
    @. u = analytical_solution.(0)
    # evolution
    for i = 1:4
        Infino.ODESolver.rk4!(example_ode!, gfs.levs[1])
    end
    t = g.levs[1].time
    @test isapprox(gfs.levs[1].u[1][1], analytical_solution(t); rtol = 1e-4)
end

@testset "RK4 Method (for New Subcycling) Tests" begin
    g = Infino.Basic.Grid(2, [[-1.0, 1.0]], 0, 0; cfl = 0.25, verbose = false)
    gfs = Infino.Basic.GridFunction(1, g)
    # initial data
    u = gfs.levs[1].u[1]
    x = gfs.levs[1].x
    @. u = analytical_solution.(0)
    # evolution
    for i = 1:4
        Infino.ODESolver.rk4_new!(example_ode!, gfs.levs[1])
    end
    t = g.levs[1].time
    @test isapprox(gfs.levs[1].u[1][1], analytical_solution(t); rtol = 1e-4)
end
