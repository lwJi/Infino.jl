module Derivs

function derivs_1st!(du, u, dx, ord)
  nx = length(u)
  istart = 1 + div(ord, 2)
  iend = nx - div(ord, 2)

  if (ord == 2)
    for i = istart:iend
      du[i] = (-u[i-1] + u[i+1]) / (2*dx)
    end
  elseif (ord == 4)
    for i = istart:iend
      du[i] = (u[i-2] - 8*u[i-1] + 8*u[i+1] - u[i+2]) / (12*dx)
    end
  else
    println("Finite difference order not supported yet: ord = ", ord)
    exit()
  end
end

function derivs_2nd!(ddu, u, dx, ord)
  nx = length(u)
  istart = 1 + div(ord, 2)
  iend = nx - div(ord, 2)

  if (ord == 2)
    for i = istart:iend
      ddu[i] = (u[i-1] - 2*u[i] + u[i+1]) / (dx*dx)
    end
  elseif (ord == 4)
    for i = istart:iend
      ddu[i] = (-u[i-2] + 16*u[i-1] - 30*u[i] + 16*u[i+1] - u[i+2]) / (12*dx*dx)
    end
  else
    println("Finite difference order not supported yet: ord = ", ord)
    exit()
  end
end

function calc_du!(du::Array{Array{Float64,1},1}, u::Array{Array{Float64,1},1},
    dx, ord=4)
  for i = 1:length(u)
    derivs_1st!(du[i], u[i], dx, ord)
  end
end

function calc_ddu!(ddu::Array{Array{Float64,1},1}, u::Array{Array{Float64,1},1},
    dx, ord=4)
  for i = 1:length(u)
    derivs_2nd!(ddu[i], u[i], dx, ord)
  end
end

end # module Derivs
