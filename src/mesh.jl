#=
"""
    mesh(start, stop, n)

Uniform mesh.

"""
function mesh(start, stop, length)
    coor = range.(start, stop, length)
    eye = Ones{Bool}.(length)
    _mesh(coor, eye)
end

_mesh(coor::NTuple{1}, eye::NTuple{1}) = coor

_mesh(coor::NTuple{2}, eye::NTuple{2}) =
    ApplyArray(kron, eye[2], coor[1]),
    ApplyArray(kron, coor[2], eye[1])
          
_mesh(coor::NTuple{3}, eye::NTuple{3}) =
    ApplyArray(kron, eye[3], eye[2], coor[1]),
    ApplyArray(kron, eye[3], coor[2], eye[1]),
    ApplyArray(kron, coor[3], eye[2], eye[1])
=#

"""
"""
mesh
#
#mesh(n::NTuple{N,Int}, args...) where {N} =
#    mesh(ntuple(i -> dirichlet, Val(N)), n, args...)

function mesh(bc, st, n, args...)
    coor = _mesh.(bc, st, n, args...)
    n̅ = @. n + 2width(st) * has(bc)
    eye = Ones{Bool}.(n̅)
    x̅, x̅ₘ = _mesh(coor, eye)
    n̅, x̅, x̅ₘ
end

function _mesh(bc::Dirichlet, st, n::Int, func=identity)
    w = width(st)
    func.(range(-w, n+w-1, n+2w) ./ (n-1)),
    func.(range(1-2w, 2(n+w)-1, n+2w) ./ (2n-2))
end

_mesh(::Periodic, ::Any, n::Int, func=identity) =
    func.(range(0, n-1, n) ./ n),
    func.(range(1, 2n-1, n) ./ 2n)

_mesh(coor::NTuple{1}, eye::NTuple{1}) =
    (first(coor[1]),),
    (last(coor[1]),)

_mesh(coor::NTuple{2}, eye::NTuple{2}) =
    (ApplyArray(kron, eye[2], first(coor[1])),
     ApplyArray(kron, first(coor[2]), eye[1])),
    (ApplyArray(kron, eye[2], last(coor[1])),
     ApplyArray(kron, last(coor[2]), eye[1]))

_mesh(coor::NTuple{3}, eye::NTuple{3}) =
    (ApplyArray(kron, eye[3], eye[2], first(coor[1])),
     ApplyArray(kron, eye[3], first(coor[2]), eye[1]),
     ApplyArray(kron, first(coor[3]), eye[2], eye[1])),
    (ApplyArray(kron, eye[3], eye[2], last(coor[1])),
     ApplyArray(kron, eye[3], last(coor[2]), eye[1]),
     ApplyArray(kron, last(coor[3]), eye[2], eye[1]))

