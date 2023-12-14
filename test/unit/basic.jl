using Infino

@testset "GridFunction" begin
    nx = 101
    ngh = 2
    nbuf = ngh * 4
    bnds = [[-4.0, 4.0], [0.6, 1.4], [0.8, 1.2]]
    g = Infino.Basic.Grid(nx, bnds, ngh, nbuf; verbose = false)
    @test isapprox(
        [g.levs[1].dx, g.levs[2].dx, g.levs[3].dx],
        [0.08, 0.04, 0.02];
        atol = 1e-12,
    )
    gf = Infino.Basic.GridFunction(1, g)
    @test isapprox(
        [gf.levs[1].x[1+nbuf], gf.levs[1].x[g.levs[1].nxa-nbuf]],
        [-4.0, 4.0];
        atol = 1e-12,
    )
    @test isapprox(
        [gf.levs[2].x[1+nbuf], gf.levs[2].x[g.levs[2].nxa-nbuf]],
        [0.6, 1.4];
        atol = 1e-12,
    )
    @test isapprox(
        [gf.levs[3].x[1+nbuf], gf.levs[3].x[g.levs[3].nxa-nbuf]],
        [0.8, 1.2];
        atol = 1e-12,
    )
end
