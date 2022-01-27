using LinearAlgebra
using SparseArrays
using CartesianCutCell

bc = (dirichlet, periodic)
st = stencil(bc)
n = (3, 4)

n̅, x̅, x̅ₘ = mesh(bc, st, n)

δ̅₋, δ̅₊ = backwarddiff(bc, n̅), forwarddiff(bc, n̅)
σ̅₋, σ̅₊ = backwardinterp(bc, n̅), forwardinterp(bc, n̅)

Δ̅ = -mapreduce(*, +, δ̅₋, δ̅₊)
Δ = restrict(Δ̅, n̅, n) |> sparse

x = restrict.(x̅, Ref(n̅), Ref(n))
xₘ = restrict.(x̅ₘ, Ref(n̅), Ref(n))

f(u, v) = π * u ^ 2 * v - √2 * v ^ 2 - √3u + 1

y = f.(x...) |> collect
b = Δ * y
yₕ = Δ \ b

@show norm(y - yₕ)

#=
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
=#
