using SparseArrays
using CartesianCutCell

#=
const center = Val(0x0)
const face = (Val(0x1), Val(0x2), Val(0x4))
const edge = (Val(0x6), Val(0x5), Val(0x3))
const node = Val(0x7)
=#

n = (3, 3)
bc = (dirichlet, periodic)

x̅, x̅ₘ = mesh(bc, n)

δ̅₋, δ̅₊ = backwarddiff(bc, n), forwarddiff(bc, n)
σ̅₋, σ̅₊ = backwardinterp(bc, n), forwardinterp(bc, n)

#μ, ρ = mask(n), pad(n)
#Δ = ρ - sparse(μ) * mapreduce(*, +, δ₋, δ₊) * μ

Δ̅ = -mapreduce(*, +, δ̅₋, δ̅₊)

Δ = prune(bc, n, Δ̅)
b = rand(prod(n))

y = Δ \ b

#

x = prune.(Ref(bc), Ref(n), x̅)
xₘ = prune.(Ref(bc), Ref(n), x̅ₘ)

using LinearAlgebra

circle(x...; x₀=zero.(x), r=one(eltype(x))/2, f=identity) =
    f(norm(r)) - f(norm(x .- x₀))

g̅ = map(circle, x̅...)

xs = map(x̅) do el
    reshape(el, @.(n+2length(bc))...)
end

#=
using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1])
surface!(ax,
hm = heatmap!(ax, randn(20, 20))
Colorbar(fig[1, 2], hm)
save("coucou.svg", fig)
=#

#map(circle, x...)
#using FillArrays
#Ones{Float64}.(n)
