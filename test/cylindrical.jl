using LinearAlgebra
using SparseArrays
using CartesianCutCell
using Vofinit

bc = (periodic, periodic, periodic)
st = stencil(bc)
n = (64, 64, 64)

n̅, x̅, x̅ₘ, dx̅ = mesh(bc, st, n)

δ̅₋, δ̅₊ = backwarddiff(bc, n̅), forwarddiff(bc, n̅)
σ̅₋, σ̅₊ = backwardinterp(bc, n̅), forwardinterp(bc, n̅)

circle(x...; x₀=zero.(x), r=one(eltype(x))/2, f=identity) =
    f(norm(r)) - f(norm(x .- x₀))

V̅ = map(tuple.(x̅...), tuple.(dx̅...)) do x₀, h₀
    getcelltype(circle, x₀, h₀)
end

V = restrict(V̅, n̅, n)

# -1 : =0 (interface)
# 0 : >0 (inside)
# 1 : <0 (outside)
cnt = count.((==).((-1, 0, 1)), Ref(V))
cnt[2] / cnt[3] * 48 / π

#V̅, A̅ = integrate(circle, x̅, dx̅)

