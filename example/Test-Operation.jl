struct LevelFunction
    nd::Int64
    nx::Int64
    u::Array{Array{Float64,1},1}
    u_p::Array{Array{Float64,1},1}

    function LevelFunction(nd, nx)
        u = Array{Array{Float64,1},1}(undef, nd)
        u_p = Array{Array{Float64,1},1}(undef, nd)
        for i = 1:nd
            u[i] = zeros(Float64, nx)
            u_p[i] = zeros(Float64, nx)
        end
        new(nd, nx, u, u_p)
    end
end

function main()
    lf = LevelFunction(2, 5)
    u = lf.u
    u_p = lf.u_p
    psi = lf.u[1]
    Pi = lf.u[2]
    @. psi = 1.0
    @. Pi = 2.0
    println("u   = ", lf.u)
    println("u_p = ", lf.u_p)
    println()
    @. u_p += u * 10.0
    println("u   = ", lf.u)
    println("u_p = ", lf.u_p)
    println()
    @. u_p = u
    println("u   = ", lf.u)
    println("u_p = ", lf.u_p)
    println()
    @. u += u_p * 100
    println("u   = ", lf.u)
    println("u_p = ", lf.u_p)
end

main()
