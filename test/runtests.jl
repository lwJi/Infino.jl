using Test, SafeTestsets

@time begin
  @safetestset "scalarwave" begin
    include("../src/example/Test.jl")
  end
end
