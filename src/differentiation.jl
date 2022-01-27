"""
"""
forwarddiff

function forwarddiff(bc, n)
    opn = _forwarddiff.(bc, n)
    eye = Eye{Bool}.(n)
    _forwarddiff(opn, eye)
end

_forwarddiff(::Dirichlet, n::Int) =
    Bidiagonal(-ones(n), ones(n-1), :U)

_forwarddiff(::Periodic, n::Int) =
    spdiagm(1-n => ones(1),
            0 => -ones(n),
            1 => ones(n-1))

_forwarddiff(opn::NTuple{1}, eye::NTuple{1}) = opn

_forwarddiff(opn::NTuple{2,Any}, eye::NTuple{2,Any}) =
    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])

_forwarddiff(opn::NTuple{3,Any}, eye::NTuple{3,Any}) =
    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])

"""
"""
backwarddiff

function backwarddiff(bc, n)
    opn = _backwarddiff.(bc, n)
    eye = Eye{Bool}.(n)
    _backwarddiff(opn, eye)
end

_backwarddiff(::Dirichlet, n::Int) =
    Bidiagonal(ones(n), -ones(n-1), :L)

_backwarddiff(::Periodic, n::Int) =
    spdiagm(-1 => -ones(n-1),
            0 => ones(n),
            n - 1 => -ones(1))

_backwarddiff(opn::NTuple{1}, eye::NTuple{1}) = opn

_backwarddiff(opn::NTuple{2,Any}, eye::NTuple{2,Any}) =
    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])

_backwarddiff(opn::NTuple{3,Any}, eye::NTuple{3,Any}) =
    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])

