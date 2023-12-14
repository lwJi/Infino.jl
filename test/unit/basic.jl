using Infino

@testset "GridFunction" begin
    tol = 1e-12
    nx = 101
    ngh = 2
    nbuf = ngh * 4
    bnds = [[-4.0, 4.0], [0.6, 1.4], [0.8, 1.2]]
    g = Infino.Basic.Grid(nx, bnds, ngh, nbuf; verbose = false)
    @test abs(g.levs[1].dx - 0.08) + abs(g.levs[2].dx - 0.04) + abs(g.levs[3].dx - 0.02) <
          tol

    gf = Infino.Basic.GridFunction(1, g)
    @test abs(gf.levs[1].x[1+nbuf] + 4.0) + abs(gf.levs[1].x[g.levs[1].nxa-nbuf] - 4.0) <
          tol
    @test abs(gf.levs[2].x[1+nbuf] - 0.6) + abs(gf.levs[2].x[g.levs[2].nxa-nbuf] - 1.4) <
          tol
    @test abs(gf.levs[3].x[1+nbuf] - 0.8) + abs(gf.levs[3].x[g.levs[3].nxa-nbuf] - 1.2) <
          tol
end
