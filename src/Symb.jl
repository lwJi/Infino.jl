module Symb

using Symbolics

#===============================================================================
RK's interpolant, termed dense output
===============================================================================#
@variables theta, h, k[1:4], yn

b = [
    theta - 3 // 2 * theta^2 + 2 // 3 * theta^3,
    theta^2 - 2 // 3 * theta^3,
    theta^2 - 2 // 3 * theta^3,
    -1 // 2 * theta^2 + 2 // 3 * theta^3,
]

y_symb = yn + sum([b[i] * k[i] for i = 1:length(b)])

dy_symb(ord) = expand_derivatives((Differential(theta)^ord)(y_symb)) / h^ord

# build function
y = build_function(y_symb, theta, yn, k, expression = Val{false})
dy(ord) = build_function(dy_symb(ord), theta, h, k, expression = Val{false})

end
