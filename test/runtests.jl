using Infino
using Test
using SafeTestsets

@time begin
    @safetestset "Basic" begin
        include("unit/basic.jl")
        include("unit/odesolver.jl")
    end
    @safetestset "Derivs" begin
        include("unit/derivs.jl")
    end
    @safetestset "ODESolver" begin
        include("unit/odesolver.jl")
    end
end
