function prune(bc, n, x::AbstractVector)
    m = length.(bc)
    cart = CartesianIndices(@. UnitRange(1+m, n+m))
    lin = LinearIndices(@. n + 2m)
    ind = reshape(lin[cart], prod(n))
    collect(x[ind])
end

function prune(bc, n, A::AbstractMatrix)
    m = length.(bc)
    cart = CartesianIndices(@. UnitRange(1+m, n+m))
    lin = LinearIndices(@. n + 2m)
    ind = reshape(lin[cart], prod(n))
    @show ind
    sparse(A[ind, ind])
end

