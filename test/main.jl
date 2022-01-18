using CartesianCutCell
using SparseArrays

#=
const center = Val(0x0)
const face = (Val(0x1), Val(0x2), Val(0x4))
const edge = (Val(0x6), Val(0x5), Val(0x3))
const node = Val(0x7)
=#

n = (3, 3, 3)
#
#a, b = (0, 0),  (1, 1)
#x = mesh(a, b, n)

δ₋, δ₊ = backwarddiff(n), forwarddiff(n)
σ₋, σ₊ = backwardinterp(n), forwardinterp(n)

μ, ρ = mask(n), pad(n)

Δ = ρ - sparse(μ) * mapreduce(*, +, δ₋, δ₊) * μ

A = sparse(Δ)
b = rand(prod(n))

x = A \ b

#using FillArrays
#Ones{Float64}.(n)
