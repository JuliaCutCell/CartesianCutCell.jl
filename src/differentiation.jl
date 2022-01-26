"""
"""
forwarddiff
#
#forwarddiff(n::NTuple{N}) where {N} =
#    forwarddiff(ntuple(i -> periodic, Val(N)), n)
#
#forwarddiff((bc,), (n,)) = (forwarddiff(bc, n),)

function forwarddiff(bc, n)
    opn = _forwarddiff.(bc, n)
    eye = @. Eye{Bool}(n + 2length(bc))
    _forwarddiff(opn, eye)
end

function _forwarddiff(bc::Dirichlet, n::Int)
    m = length(bc)
    Bidiagonal(-ones(n+2m), ones(n+2m-1), :U)
end

function _forwarddiff(bc::Periodic, n::Int)
    m = length(bc)
    spdiagm(1-n => ones(1),
            0 => -ones(n+2m),
            1 => ones(n+2m-1))
end

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
#
#backwarddiff(n::NTuple{N}) where {N} =
#    backwarddiff(ntuple(i -> periodic, Val(N)), n)
#
#backwarddiff((bc,), (n,)) = (backwarddiff(bc, n),)

function backwarddiff(bc, n)
    opn = _backwarddiff.(bc, n)
    eye = @. Eye{Bool}(n + 2length(bc))
    _backwarddiff(opn, eye)
end

function _backwarddiff(bc::Dirichlet, n::Int)
    m = length(bc)
    Bidiagonal(ones(n+2m), -ones(n+2m-1), :L)
end

function _backwarddiff(bc::Periodic, n::Int)
    m = length(bc)
    spdiagm(-1 => -ones(n+2m-1),
            0 => ones(n+2m),
            n - 1 => -ones(1))
end

_backwarddiff(opn::NTuple{1}, eye::NTuple{1}) = opn

_backwarddiff(opn::NTuple{2,Any}, eye::NTuple{2,Any}) =
    ApplyArray(kron, eye[2], opn[1]),
    ApplyArray(kron, opn[2], eye[1])

_backwarddiff(opn::NTuple{3,Any}, eye::NTuple{3,Any}) =
    ApplyArray(kron, eye[3], eye[2], opn[1]),
    ApplyArray(kron, eye[3], opn[2], eye[1]),
    ApplyArray(kron, opn[3], eye[2], eye[1])

