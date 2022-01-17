"""
"""
forwarddiff

forwarddiff(n::NTuple{N}) where {N} =
    forwarddiff(ntuple(i -> dirichlet, Val(N)), n)

forwarddiff((bc,), (n,)) = (forwarddiff(bc, n),)

forwarddiff(::Dirichlet, n::Int) =
    Bidiagonal(-ones(n), ones(n-1), :U)

function forwarddiff(::Periodic, n::Int)
    spdiagm(1-n => ones(1),
            0 => -ones(n),
            1 => ones(n-1))
end

function forwarddiff(bc::NTuple{2}, n::NTuple{2})
    opn = forwarddiff.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])
end

function forwarddiff(bc::NTuple{3}, n::NTuple{3})
    opn = forwarddiff.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])
end

"""
"""
backwarddiff

backwarddiff(n::NTuple{N}) where {N} =
    backwarddiff(ntuple(i -> dirichlet, Val(N)), n)

backwarddiff((bc,), (n,)) = (backwarddiff(bc, n),)

backwarddiff(::Dirichlet, n::Int) =
    Bidiagonal(ones(n), -ones(n-1), :L)

function backwarddiff(::Periodic, n::Int)
    spdiagm(-1 => -ones(n-1),
            0 => ones(n),
            n - 1 => -ones(1))
end

function backwarddiff(bc::NTuple{2}, n::NTuple{2})
    opn = backwarddiff.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])
end

function backwarddiff(bc::NTuple{3}, n::NTuple{3})
    opn = backwarddiff.(bc, n)
    eye = I.(n)

    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])
end

