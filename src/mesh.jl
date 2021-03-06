"""
"""
mesh

function mesh(bc, st, n, args...)
    coor = _mesh.(bc, st, n, args...)
    n̅ = @. n + 2width(st) * has(bc)
    eye = Ones{Bool}.(n̅)
    x̅, x̅ₘ, dx̅ = _mesh(coor, eye)
    n̅, x̅, x̅ₘ, dx̅
end

function _mesh(bc::Dirichlet, st, n::Int, func=identity)
    w = width(st)
    func.(range(-w, n+w-1, n+2w) ./ (n-1)),
    func.(range(1-2w, 2(n+w)-1, n+2w) ./ (2n-2)),
    func.(range(1-w, n+w, n+2w) ./ (n-1)) .-
    func.(range(-w, n+w-1, n+2w) ./ (n-1))
end

_mesh(::Periodic, ::Any, n::Int, func=identity) =
    func.(range(0, n-1, n) ./ n),
    func.(range(1, 2n-1, n) ./ 2n),
    func.(range(1, n, n) ./ n) .-
    func.(range(0, n-1, n) ./ n)

_mesh(coor::NTuple{1}, eye::NTuple{1}) =
    (coor[1][1],),
    (coor[1][2],),
    (coor[1][3],)

_mesh(coor::NTuple{2}, eye::NTuple{2}) =
    (ApplyArray(kron, eye[2], coor[1][1]),
     ApplyArray(kron, coor[2][1], eye[1])),
    (ApplyArray(kron, eye[2], coor[1][2]),
     ApplyArray(kron, coor[2][2], eye[1])),
    (ApplyArray(kron, eye[2], coor[1][3]),
     ApplyArray(kron, coor[2][3], eye[1]))

_mesh(coor::NTuple{3}, eye::NTuple{3}) =
    (ApplyArray(kron, eye[3], eye[2], coor[1][1]),
     ApplyArray(kron, eye[3], coor[2][1], eye[1]),
     ApplyArray(kron, coor[3][1], eye[2], eye[1])),
    (ApplyArray(kron, eye[3], eye[2], coor[1][2]),
     ApplyArray(kron, eye[3], coor[2][2], eye[1]),
     ApplyArray(kron, coor[3][2], eye[2], eye[1])),
    (ApplyArray(kron, eye[3], eye[2], coor[1][3]),
     ApplyArray(kron, eye[3], coor[2][3], eye[1]),
     ApplyArray(kron, coor[3][3], eye[2], eye[1]))

