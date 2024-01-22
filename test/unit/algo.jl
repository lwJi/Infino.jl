using Infino

@testset "1st-Order Interpolaton" begin
    ord = 1
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-2)
end

@testset "2nd-Order Interpolaton" begin
    ord = 2
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-4)
end

@testset "3rd-Order Interpolaton" begin
    ord = 3
    x = LinRange(0.1, 0.4, 4)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 2, ord), exp(0.25); rtol = 1e-5)
end

@testset "5th-Order Interpolaton" begin
    ord = 5
    x = LinRange(0.1, 0.4, 6)
    y = exp.(x)
    @test isapprox(Infino.Algo.Interpolation(y, 3, ord), exp(0.25); rtol = 1e-6)
end
