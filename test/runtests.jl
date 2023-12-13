using Infino
using Test
using SafeTestsets

@time begin
  @safetestset "Basic" begin
    include("unit/basis.jl")
  end
end
