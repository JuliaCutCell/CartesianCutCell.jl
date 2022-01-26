using CartesianCutCell
using SparseArrays

#=
const center = Val(0x0)
const face = (Val(0x1), Val(0x2), Val(0x4))
const edge = (Val(0x6), Val(0x5), Val(0x3))
const node = Val(0x7)
=#

n = (4, 3)
bc = (dirichlet, periodic)

x, xm = mesh(bc, n)

δ₋, δ₊ = backwarddiff(bc, n), forwarddiff(bc, n)
σ₋, σ₊ = backwardinterp(bc, n), forwardinterp(bc, n)
#
#μ, ρ = mask(n), pad(n)

#Δ = ρ - sparse(μ) * mapreduce(*, +, δ₋, δ₊) * μ
Δ = -mapreduce(*, +, δ₋, δ₊)

v = rand(prod(@. n + 2length(bc)))
prune(bc, n, v)

A = prune(bc, n, Δ)
b = rand(prod(n))

x = A \ b

#using FillArrays
#Ones{Float64}.(n)
