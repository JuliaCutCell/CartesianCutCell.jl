"""
"""
mask

mask(n::NTuple{N}) where {N} =
    mask(ntuple(i -> dirichlet, Val(N)), n)

function mask(bc::NTuple{N}, n::NTuple{N}) where {N}
    opn = _mask.(bc, n)

    ntuple(Val(N)) do i

    end
    center = kron(first.(opn)...)

    center, face, edge, node =
end

function mask(bc::NTuple{1}, n::NTuple{1})
    opn = _mask.(bc, n)
    opn[1]
end

function mask(bc::NTuple{2}, n::NTuple{2})
    opn = _mask.(bc, n)
    ApplyArray(kron, opn[2], opn[1])
end

function mask(bc::NTuple{3}, n::NTuple{3})
    opn = _mask.(bc, n)
    ApplyArray(kron, opn[3], opn[2], opn[1])
end

_mask(::Dirichlet, n::Int) =
    (Diagonal(BitVector(i != n for i in 1:n)),
     Eye{Bool}(n))

_mask(::Periodic, n::Int) =
    (Eye{Bool}(n),
     Eye{Bool}(n))

