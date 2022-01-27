function restrict(x̅::AbstractVector, n̅, n)
    s = @. (n̅ - n) ÷ 2
    cart = CartesianIndices(@. UnitRange(1+s, n+s))
    lin = LinearIndices(@. n + 2s)
    ind = reshape(lin[cart], prod(n))
    view(x̅, ind)
end

function restrict(A̅::AbstractMatrix, n̅, n)
    s = @. (n̅ - n) ÷ 2
    cart = CartesianIndices(@. UnitRange(1+s, n+s))
    lin = LinearIndices(@. n + 2s)
    ind = reshape(lin[cart], prod(n))
    view(A̅, ind, ind)
end

