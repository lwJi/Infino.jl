module InitialData

include("Derivs.jl")

function Gaussian!(gfs)
    amp = 1.0
    sig = 0.2
    x0 = 0.0

    for l = 1:length(gfs.levs)
        psi = gfs.levs[l].u[1]
        Pi = gfs.levs[l].u[2]
        x = gfs.levs[l].x

        @. psi = amp * exp(-((x - x0) / sig)^2)
        @. Pi = 0.0
    end
end

end # module InitialData
