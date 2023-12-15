using Infino
using Test
using SafeTestsets

@time begin
    @safetestset "Algo" begin
        include("unit/algo.jl")
    end
    @safetestset "Basic" begin
        include("unit/basic.jl")
    end
    @safetestset "Derivs" begin
        include("unit/derivs.jl")
    end
    @safetestset "ODESolver" begin
        include("unit/odesolver.jl")
    end
    @safetestset "ScalarWave" begin
        include("unit/scalarwave.jl")
    end
end
