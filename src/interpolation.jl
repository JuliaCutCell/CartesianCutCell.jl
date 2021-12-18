"""
"""
forwardinterp

forwardinterp(n::NTuple{N}) where {N} =
    forwardinterp(ntuple(i -> dirichlet, Val(N)), n)

forwardinterp((bc,), (n,)) = (forwardinterp(bc, n),)

function forwardinterp(::Dirichlet, n::Int)
    spdiagm(0 => ones(n) / 2,
            1 => ones(n-1) / 2)
end

function forwardinterp(::Periodic, n::Int)
    spdiagm(1-n => ones(1) / 2,
            0 => ones(n) / 2,
            1 => ones(n-1) / 2)
end

function forwardinterp(bc::NTuple{2}, n::NTuple{2})
    opn = forwardinterp.(bc, n)
    eye = I.(n)

    kron(eye[2], opn[1]),
    kron(opn[2], eye[1])
end

function forwardinterp(bc::NTuple{3}, n::NTuple{3})
    opn = forwardinterp.(bc, n)
    eye = I.(n)

    kron(eye[3], eye[2], opn[1]),
    kron(eye[3], opn[2], eye[1]),
    kron(opn[3], eye[2], eye[1])
end

"""
"""
backwardinterp

backwardinterp(n::NTuple{N}) where {N} =
    backwardinterp(ntuple(i -> dirichlet, Val(N)), n)

backwardinterp((bc,), (n,)) = (backwardinterp(bc, n),)

function backwardinterp(::Dirichlet, n::Int)
    spdiagm(1 => ones(n-1) / 2,
            0 => ones(n) / 2)
end

function backwardinterp(::Periodic, n::Int)
    spdiagm(-1 => ones(n-1) / 2,
            0 => ones(n) / 2,
            n - 1 => ones(1) / 2)
end

function backwardinterp(bc::NTuple{2}, n::NTuple{2})
    opn = backwardinterp.(bc, n)
    eye = I.(n)

    kron(eye[2], opn[1]),
    kron(opn[2], eye[1])
end

function backwardinterp(bc::NTuple{3}, n::NTuple{3})
    opn = backwardinterp.(bc, n)
    eye = I.(n)

    kron(eye[3], eye[2], opn[1]),
    kron(eye[3], opn[2], eye[1]),
    kron(opn[3], eye[2], eye[1])
end

