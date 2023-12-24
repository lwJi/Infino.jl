include("Symb.jl")
using .Symb
using JuliaFormatter
using Symbolics

f_y = build_function(y, theta, yn, k, target = Symbolics.JuliaTarget())
f_dy(ord) = build_function(dy(ord), theta, h, k, target = Symbolics.JuliaTarget())

fname = "../src/DenseOutput.jl"
open(fname, "w") do file
    println(file, "module DenseOutput")
    println(file, "")
    println(file, "y = ")
    println(file, string(f_y))
    println(file, "")
    for i = 1:3
        println(file, "dy$i = ")
        println(file, string(f_dy(i)))
        println(file, "")
    end
    println(file, "end")
end
format(fname)
