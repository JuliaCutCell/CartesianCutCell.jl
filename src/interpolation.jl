"""
"""
forwardinterp
#
#forwardinterp(n::NTuple{N}) where {N} =
#    forwardinterp(ntuple(i -> periodic, Val(N)), n)
#
#forwardinterp((bc,), (n,)) = (forwardinterp(bc, n),)

function forwardinterp(bc, n)
    opn = _forwardinterp.(bc, n)
    eye = Eye{Bool}.(n)
    _forwardinterp(opn, eye)
end

_forwardinterp(::Dirichlet, n::Int) =
    Bidiagonal(ones(n) / 2, ones(n-1) / 2, :U)

_forwardinterp(::Periodic, n::Int) =
    spdiagm(1-n => ones(1) / 2,
            0 => ones(n) / 2,
            1 => ones(n-1) / 2)

_forwardinterp(opn::NTuple{1}, eye::NTuple{1}) = opn

_forwardinterp(opn::NTuple{2,Any}, eye::NTuple{2,Any}) =
    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])

_forwardinterp(opn::NTuple{3,Any}, eye::NTuple{3,Any}) =
    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])

"""
"""
backwardinterp
#
#backwardinterp(n::NTuple{N}) where {N} =
#    backwardinterp(ntuple(i -> periodic, Val(N)), n)
#
#backwardinterp((bc,), (n,)) = (backwardinterp(bc, n),)

function backwardinterp(bc, n)
    opn = _backwardinterp.(bc, n)
    eye = Eye{Bool}.(n)
    _backwardinterp(opn, eye)
end

_backwardinterp(::Dirichlet, n::Int) =
    Bidiagonal(ones(n) / 2, ones(n-1) / 2, :L)

_backwardinterp(::Periodic, n::Int) =
    spdiagm(-1 => ones(n-1) / 2,
            0 => ones(n) / 2,
            n - 1 => ones(1) / 2)

_backwardinterp(opn::NTuple{1}, eye::NTuple{1}) = opn

_backwardinterp(opn::NTuple{2,Any}, eye::NTuple{2,Any}) =
    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])

_backwardinterp(opn::NTuple{3,Any}, eye::NTuple{3,Any}) =
    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])

