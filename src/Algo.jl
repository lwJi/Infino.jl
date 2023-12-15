module Algo

function Interpolation(u, i, ord)
    if ord == 1
        return (u[i] + u[i+1]) * 0.5
    elseif ord == 2
        return (-u[i-1] + 6 * u[i] + 3 * u[i+1]) * 0.125
    elseif ord == 3
        return (-u[i-1] + 9 * u[i] + 9 * u[i+1] - u[i+2]) * 0.0625
    else
        println("Interpolation order not supported yet: ord = ", ord)
        exit()
    end
end

end
