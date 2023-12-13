using Infino

@testset "Derivs" begin
  tol = 1e-8
  nx = 11
  x = LinRange(-0.5, 0.5, nx)
  dx = x[2] - x[1]
  y1 = zeros(Float64, nx)
  y2 = zeros(Float64, nx)
  dy1 = zeros(Float64, nx)
  ddy2 = zeros(Float64, nx)
  @. y1 = 2.0 * x
  @. y2 = x * x
  Infino.Derivs.derivs_1st!(dy1, y1, dx, 2)
  @test all(dy1[3:nx-2] .- 2.0 .< tol)
  Infino.Derivs.derivs_1st!(dy1, y1, dx, 4)
  @test all(dy1[3:nx-2] .- 2.0 .< tol)
  Infino.Derivs.derivs_2nd!(ddy2, y2, dx, 2)
  @test all(ddy2[3:nx-2] .- 2.0 .< tol)
  Infino.Derivs.derivs_2nd!(ddy2, y2, dx, 4)
  @test all(ddy2[3:nx-2] .- 2.0 .< tol)

end
