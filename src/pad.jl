"""
"""
pad

pad(n::NTuple{N}) where {N} =
    pad(ntuple(i -> dirichlet, Val(N)), n)

pad((bc,), (n,)) = pad(bc, n)

pad(::Dirichlet, n::Int) =
    Diagonal(BitVector(i == n for i in 1:n))

pad(::Periodic, n::Int) = Eye{Bool}(n)

function pad(bc::NTuple{2}, n::NTuple{2})
    opn = pad.(bc, n)
    eye = Eye{Bool}.(n)
    ApplyArray(kron, eye[2], opn[1]) .|
    ApplyArray(kron, opn[2], eye[1])
end

function pad(bc::NTuple{3}, n::NTuple{3})
    opn = pad.(bc, n)
    eye = Eye{Bool}.(n)
    ApplyArray(kron, eye[3], eye[2], opn[1]) .|
    ApplyArray(kron, eye[3], opn[2], eye[1]) .|
    ApplyArray(kron, opn[3], eye[2], eye[1])
end

