"""
"""
forwardinterp

forwardinterp(n::NTuple{N}) where {N} =
    forwardinterp(ntuple(i -> dirichlet, Val(N)), n)

forwardinterp((bc,), (n,)) = (forwardinterp(bc, n),)

forwardinterp(::Dirichlet, n::Int) =
    Bidiagonal(ones(n) / 2, ones(n-1) / 2, :U)

function forwardinterp(::Periodic, n::Int)
    spdiagm(1-n => ones(1) / 2,
            0 => ones(n) / 2,
            1 => ones(n-1) / 2)
end

function forwardinterp(bc::NTuple{2}, n::NTuple{2})
    opn = forwardinterp.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])
end

function forwardinterp(bc::NTuple{3}, n::NTuple{3})
    opn = forwardinterp.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])
end

"""
"""
backwardinterp

backwardinterp(n::NTuple{N}) where {N} =
    backwardinterp(ntuple(i -> dirichlet, Val(N)), n)

backwardinterp((bc,), (n,)) = (backwardinterp(bc, n),)

backwardinterp(::Dirichlet, n::Int) =
    Bidiagonal(ones(n) / 2, ones(n-1) / 2, :L)

function backwardinterp(::Periodic, n::Int)
    spdiagm(-1 => ones(n-1) / 2,
            0 => ones(n) / 2,
            n - 1 => ones(1) / 2)
end

function backwardinterp(bc::NTuple{2}, n::NTuple{2})
    opn = backwardinterp.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])
end

function backwardinterp(bc::NTuple{3}, n::NTuple{3})
    opn = backwardinterp.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])
end

