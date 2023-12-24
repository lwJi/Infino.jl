module DenseOutput

y = 
function (theta, yn, k)
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:373 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:374 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:375 =#
    (+)((+)((+)((+)(yn, (*)((getindex)(k, 1), (+)((+)(theta, (*)(-3//2, (^)(theta, 2))), (*)(2//3, (^)(theta, 3))))), (*)((getindex)(k, 2), (+)((^)(theta, 2), (*)(-2//3, (^)(theta, 3))))), (*)((getindex)(k, 3), (+)((^)(theta, 2), (*)(-2//3, (^)(theta, 3))))), (*)((getindex)(k, 4), (+)((*)(-1//2, (^)(theta, 2)), (*)(2//3, (^)(theta, 3)))))
end

dy1 = 
function (theta, h, k)
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:373 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:374 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:375 =#
    (/)((+)((+)((+)((*)((getindex)(k, 1), (+)((+)(1, (*)(-3//1, theta)), (*)(2//1, (^)(theta, 2)))), (*)((getindex)(k, 2), (+)((*)(2, theta), (*)(-2//1, (^)(theta, 2))))), (*)((getindex)(k, 3), (+)((*)(2, theta), (*)(-2//1, (^)(theta, 2))))), (*)((getindex)(k, 4), (+)((*)(-1//1, theta), (*)(2//1, (^)(theta, 2))))), h)
end

dy2 = 
function (theta, h, k)
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:373 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:374 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:375 =#
    (/)((+)((+)((+)((*)((getindex)(k, 1), (+)(-3//1, (*)(4//1, theta))), (*)((getindex)(k, 2), (+)(2, (*)(-4//1, theta)))), (*)((getindex)(k, 3), (+)(2, (*)(-4//1, theta)))), (*)((getindex)(k, 4), (+)(-1//1, (*)(4//1, theta)))), (^)(h, 2))
end

dy3 = 
function (theta, h, k)
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:373 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:374 =#
    #= /Users/liwei/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:375 =#
    (/)((+)((+)((+)((*)(4//1, (getindex)(k, 1)), (*)(-4//1, (getindex)(k, 2))), (*)(-4//1, (getindex)(k, 3))), (*)(4//1, (getindex)(k, 4))), (^)(h, 3))
end

end
