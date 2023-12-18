using Infino

@testset "2nd-Order Accurate 1st-Derivative" begin
    ord = 2
    ngh = Int(ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    dy = zeros(Float64, nx)
    @. y = 2.0 * x
    Infino.Derivs.derivs_1st!(dy, y, dx, ord)
    @test isapprox(dy[1+ngh:nx-ngh], ones(Float64, nx - 2 * ngh) * 2.0; rtol = 1e-12)
end

@testset "4th-Order Accurate 1st-Derivative" begin
    ord = 4
    ngh = Int(ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    dy = zeros(Float64, nx)
    @. y = 2.0 * x
    Infino.Derivs.derivs_1st!(dy, y, dx, ord)
    @test isapprox(dy[1+ngh:nx-ngh], ones(Float64, nx - 2 * ngh) * 2.0; rtol = 1e-12)
end

@testset "2nd-Order Accurate 2nd-Derivative" begin
    ord = 2
    ngh = Int(ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    ddy = zeros(Float64, nx)
    @. y = x * x
    Infino.Derivs.derivs_2nd!(ddy, y, dx, ord)
    @test isapprox(ddy[1+ngh:nx-ngh], ones(Float64, nx - 2 * ngh) * 2.0; rtol = 1e-12)
end

@testset "4th-Order Accurate 2nd-Derivative" begin
    ord = 4
    ngh = Int(ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    ddy = zeros(Float64, nx)
    @. y = x * x
    Infino.Derivs.derivs_2nd!(ddy, y, dx, ord)
    @test isapprox(ddy[1+ngh:nx-ngh], ones(Float64, nx - 2 * ngh) * 2.0; rtol = 1e-12)
end

@testset "2nd-Order Derivatie's KO Dissipation" begin
    ord = 2
    diss_ord = ord + 2
    ngh = Int(diss_ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    d4y = zeros(Float64, nx)
    @. y = x * x * x * x
    Infino.Derivs.derivs_diss!(d4y, y, dx, ord)
    @test isapprox(
        d4y[1+ngh:nx-ngh],
        -ones(Float64, nx - 2 * ngh) / 2^diss_ord * dx^(diss_ord - 1) * 24;
        rtol = 1e-12,
    )
end

@testset "4th-Order Derivatie's KO Dissipation" begin
    ord = 4
    diss_ord = ord + 2
    ngh = Int(diss_ord / 2)
    nx = 11
    x = LinRange(-0.5, 0.5, nx)
    dx = x[2] - x[1]
    y = zeros(Float64, nx)
    d4y = zeros(Float64, nx)
    @. y = x * x * x * x * x * x
    Infino.Derivs.derivs_diss!(d4y, y, dx, ord)
    @test isapprox(
        d4y[1+ngh:nx-ngh],
        ones(Float64, nx - 2 * ngh) / 2^diss_ord * dx^(diss_ord - 1) * 720;
        rtol = 1e-12,
    )
end
