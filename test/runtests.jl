using Infino
using Test
using SafeTestsets

@time begin
  @safetestset "Basic" begin
    include("unit/basic.jl")
    include("unit/derivs.jl")
  end
end
