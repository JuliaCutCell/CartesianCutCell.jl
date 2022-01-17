using CartesianCutCell
using SparseArrays

n = (3, 5)

a, b = (0, 0),  (1, 1)
x = mesh(a, b, n)

δ₋, δ₊ = backwarddiff(n), forwarddiff(n)
σ₋, σ₊ = backwardinterp(n), forwardinterp(n)

Δ = sparse(mapreduce(*, +, δ₋, δ₊))

#using FillArrays
#Ones{Float64}.(n)
