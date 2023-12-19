module Derivs

#===============================================================================
Derivatives:
    * single point
===============================================================================#
function deriv_1st(u, i, dx, ord)
    if ord == 2
        return (-u[i-1] + u[i+1]) / (2 * dx)
    elseif ord == 4
        return (u[i-2] - 8 * u[i-1] + 8 * u[i+1] - u[i+2]) / (12 * dx)
    else
        println("Finite difference order not supported yet: ord = ", ord)
        exit()
    end
end

function deriv_2nd(u, i, dx, ord)
    if ord == 2
        return (u[i-1] - 2 * u[i] + u[i+1]) / (dx * dx)
    elseif (ord == 4)
        return (-u[i-2] + 16 * u[i-1] - 30 * u[i] + 16 * u[i+1] - u[i+2]) / (12 * dx * dx)
    else
        println("Finite difference order not supported yet: ord = ", ord)
        exit()
    end
end

function deriv_diss(u, i, dx, diss_ord)
    if diss_ord == 4
        return ((u[i+2] + u[i-2]) - 4 * (u[i+1] + u[i-1]) + 6 * u[i]) / dx
    elseif diss_ord == 6
        return (
            (u[i+3] + u[i-3]) - 6 * (u[i+2] + u[i-2]) + 15 * (u[i+1] + u[i-1]) - 20 * u[i]
        ) / dx
    else
        println("KO diss order not supported yet: ord = ", ord)
        exit()
    end
end

#===============================================================================
Derivatives:
    * single level
===============================================================================#
function derivs_1st!(du, u, dx, ord)
    istart = 1 + div(ord, 2)
    iend = length(u) - div(ord, 2)
    for i = istart:iend
        du[i] = deriv_1st(u, i, dx, ord)
    end
end

function derivs_2nd!(ddu, u, dx, ord)
    istart = 1 + div(ord, 2)
    iend = length(u) - div(ord, 2)
    for i = istart:iend
        ddu[i] = deriv_2nd(u, i, dx, ord)
    end
end

function derivs_diss!(diss, u, dx, ord)
    diss_ord = ord + 2
    sign = (mod(diss_ord, 4) == 0 ? -1 : +1)
    istart = 1 + div(diss_ord, 2)
    iend = length(u) - div(diss_ord, 2)
    for i = istart:iend
        diss[i] = deriv_diss(u, i, dx, diss_ord) * (sign / 2^(diss_ord))
    end
end

#===============================================================================
Derivatives:
    * single level
    * for varlist
===============================================================================#
function calc_du!(du::Array{Array{Float64,1},1}, u::Array{Array{Float64,1},1}, dx, ord)
    for i = 1:length(u)
        derivs_1st!(du[i], u[i], dx, ord)
    end
end

function calc_ddu!(ddu::Array{Array{Float64,1},1}, u::Array{Array{Float64,1},1}, dx, ord)
    for i = 1:length(u)
        derivs_2nd!(ddu[i], u[i], dx, ord)
    end
end

end # module Derivs
