using Infino

@testset "1st-Order Interpolaton" begin
    ord = 1
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-2)
end

@testset "3rd-Order Interpolaton" begin
    ord = 3
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-5)
end
