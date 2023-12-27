using Infino

function analytical_psi(t, x)
    return sin(2 * pi * (x - t))
end

function analytical_Pi(t, x)
    return -2 * pi * cos(2 * pi * (x - t))
end

@testset "Scalar Wave Evolution on Unit Grid" begin
    g = Infino.Basic.Grid(100, [[-0.5, 0.5 - 0.01]], 3, 3; cfl = 0.25, verbose = false)
    gfs = Infino.Basic.GridFunction(2, g)
    nxa = g.levs[1].nxa
    nbuf = g.levs[1].nbuf
    # initial data
    psi = gfs.levs[1].u[1]
    Pi = gfs.levs[1].u[2]
    x = gfs.levs[1].x
    @. psi = analytical_psi(0, x)
    @. Pi = analytical_Pi(0, x)
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    # evolution
    for i = 1:4
        Infino.ODESolver.rk4!(Infino.Physical.WaveRHS!, gfs.levs[1])
    end
    t = g.levs[1].time
    @test isapprox(
        gfs.levs[1].u[1][1+nbuf:nxa-nbuf],
        analytical_psi.(t, x)[1+nbuf:nxa-nbuf];
        rtol = 1e-7,
    )
    @test isapprox(
        gfs.levs[1].u[2][1+nbuf:nxa-nbuf],
        analytical_Pi.(t, x)[1+nbuf:nxa-nbuf];
        rtol = 1e-6,
    )
end

@testset "Scalar Wave Evolution on 3 Levels Grid without Subcycling" begin
    g = Infino.Basic.Grid(
        100,
        [[-0.5, 0.5 - 0.01], [-0.25, 0.25 - 0.005], [-0.125, 0.125 - 0.0025]],
        3,
        3;
        cfl = 0.25,
        subcycling = false,
        verbose = false,
    )
    gfs = Infino.Basic.GridFunction(2, g)
    lmax = length(gfs.levs)
    # initial data
    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x
        @. psi = analytical_psi(0, x)
        @. Pi = analytical_Pi(0, x)
    end
    for l = lmax-1:-1:1
        Infino.Sync.Restriction(gfs, l)
    end
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    # evolution
    for i = 1:16
        Infino.ODESolver.Evolve!(Infino.Physical.WaveRHS!, gfs)
    end
    t = g.time
    for l = 1:lmax
        nxa = g.levs[l].nxa
        nbuf = g.levs[l].nbuf
        x = gfs.levs[l].x
        @test isapprox(
            gfs.levs[l].u[1][1+nbuf:nxa-nbuf],
            analytical_psi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-6,
        )
        @test isapprox(
            gfs.levs[l].u[2][1+nbuf:nxa-nbuf],
            analytical_Pi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-5,
        )
    end
end

@testset "Scalar Wave Evolution on 3 Levels Grid with Subcycling" begin
    g = Infino.Basic.Grid(
        100,
        [[-0.5, 0.5 - 0.01], [-0.25, 0.25 - 0.005], [-0.125, 0.125 - 0.0025]],
        3,
        3 * 4;
        cfl = 0.25,
        verbose = false,
    )
    gfs = Infino.Basic.GridFunction(2, g)
    lmax = length(gfs.levs)
    # initial data
    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x
        @. psi = analytical_psi(0, x)
        @. Pi = analytical_Pi(0, x)
    end
    for l = lmax-1:-1:1
        Infino.Sync.Restriction(gfs, l)
    end
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    Infino.InitialData.MarchBackwards!(gfs)
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    # evolution
    for i = 1:4
        Infino.ODESolver.Evolve!(Infino.Physical.WaveRHS!, gfs)
    end
    t = g.time
    for l = 1:lmax
        nxa = g.levs[l].nxa
        nbuf = g.levs[l].nbuf
        x = gfs.levs[l].x
        @test isapprox(
            gfs.levs[l].u[1][1+nbuf:nxa-nbuf],
            analytical_psi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-7,
        )
        @test isapprox(
            gfs.levs[l].u[2][1+nbuf:nxa-nbuf],
            analytical_Pi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-5,
        )
    end
end

@testset "Scalar Wave Evolution on 3 Levels Grid with Mongwane's Subcycling" begin
    g = Infino.Basic.Grid(
        100,
        [[-0.5, 0.5 - 0.01], [-0.25, 0.25], [-0.125, 0.125]],
        2,
        2;
        cfl = 0.25,
        verbose = false,
    )
    gfs = Infino.Basic.GridFunction(2, g)
    lmax = length(gfs.levs)
    # initial data
    for l = 1:lmax
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x
        @. psi = analytical_psi(0, x)
        @. Pi = analytical_Pi(0, x)
    end
    for l = lmax-1:-1:1
        Infino.Sync.Restriction(gfs, l)
    end
    Infino.Boundary.ApplyPeriodicBoundaryCondition!(gfs)
    # evolution
    for i = 1:4
        Infino.ODESolver.Evolve!(Infino.Physical.WaveRHS!, gfs; Mongwane = true)
    end
    t = g.time
    for l = 1:lmax
        nxa = g.levs[l].nxa
        nbuf = g.levs[l].nbuf
        x = gfs.levs[l].x
        @test isapprox(
            gfs.levs[l].u[1][1+nbuf:nxa-nbuf],
            analytical_psi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-5,
        )
        @test isapprox(
            gfs.levs[l].u[2][1+nbuf:nxa-nbuf],
            analytical_Pi.(t, x)[1+nbuf:nxa-nbuf];
            rtol = 1e-4,
        )
    end
end
