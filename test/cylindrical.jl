using LinearAlgebra
using SparseArrays
using CartesianCutCell

bc = (dirichlet, periodic)
st = stencil(bc)
n = (3, 4)

n̅, x̅, x̅ₘ, dx̅ = mesh(bc, st, n)

δ̅₋, δ̅₊ = backwarddiff(bc, n̅), forwarddiff(bc, n̅)
σ̅₋, σ̅₊ = backwardinterp(bc, n̅), forwardinterp(bc, n̅)

circle(x...; x₀=zero.(x), r=one(eltype(x))/2, f=identity) =
    f(norm(r)) - f(norm(x .- x₀))

V̅, A̅ = integrate(circle, x̅, dx̅)

