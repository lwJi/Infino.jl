using Infino

@testset "2nd-Order Interpolaton" begin
    ord = 2
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-2)
end

@testset "4th-Order Interpolaton" begin
    ord = 4
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-5)
end
