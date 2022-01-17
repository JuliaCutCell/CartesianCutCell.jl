"""
    mesh(start, stop, n)

Uniform mesh.

"""
function mesh(start, stop, length)
    coor = collect.(range.(start, stop, length))
    eye = ones.(Bool, length)
    eye = ones.(length)
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

