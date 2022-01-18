"""
"""
mask

mask(n::NTuple{N}) where {N} =
    mask(ntuple(i -> dirichlet, Val(N)), n)

mask((bc,), (n,)) = mask(bc, n)

mask(::Dirichlet, n::Int) =
    Diagonal(BitVector(i != n for i in 1:n))

mask(::Periodic, n::Int) = Eye{Bool}(n)

function mask(bc::NTuple{2}, n::NTuple{2})
    opn = mask.(bc, n)
    ApplyArray(kron, opn[2], opn[1])
end

function mask(bc::NTuple{3}, n::NTuple{3})
    opn = mask.(bc, n)
    ApplyArray(kron, opn[3], opn[2], opn[1])
end

