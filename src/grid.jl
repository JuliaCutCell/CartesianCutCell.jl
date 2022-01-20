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

mesh(n::NTuple{N,Int}, args...) where {N} =
    mesh(ntuple(i -> dirichlet, Val(N)), n, args...)

function mesh(bc, n, args...)
    coor = _mesh.(bc, n, args...)
    eye = Ones{Bool}.(n)
    _mesh(coor, eye)
end

_mesh(::Dirichlet, n::Int, func=identity) =
    func.(range(0, n, n) ./ n)

_mesh(::Periodic, n::Int, func=identity) =
    func.(range(0, n-1, n) ./ n)

_mesh(coor::NTuple{1}, eye::NTuple{1}) = coor

_mesh(coor::NTuple{2}, eye::NTuple{2}) =
    ApplyArray(kron, eye[2], coor[1]),
    ApplyArray(kron, coor[2], eye[1])

_mesh(coor::NTuple{3}, eye::NTuple{3}) =
    ApplyArray(kron, eye[3], eye[2], coor[1]),
    ApplyArray(kron, eye[3], coor[2], eye[1]),
    ApplyArray(kron, coor[3], eye[2], eye[1])

