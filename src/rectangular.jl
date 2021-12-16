"""
    laplacian(args...)

Operators on rectangular domains where :

1. There are no interior boundaries,
1. Exterior boundaries are treated as Dirichlet only for now.

# Examples

```jldoctest
julia> laplacian((3,))
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 7 stored entries:
  2.0  -1.0    ⋅
 -1.0   2.0  -1.0
   ⋅   -1.0   2.0

julia> laplacian((3, 2))
6×6 SparseArrays.SparseMatrixCSC{Float64, Int64} with 20 stored entries:
  4.0  -1.0    ⋅   -1.0    ⋅     ⋅
 -1.0   4.0  -1.0    ⋅   -1.0    ⋅
   ⋅   -1.0   4.0    ⋅     ⋅   -1.0
 -1.0    ⋅     ⋅    4.0  -1.0    ⋅
   ⋅   -1.0    ⋅   -1.0   4.0  -1.0
   ⋅     ⋅   -1.0    ⋅   -1.0   4.0

julia> laplacian((3, 2, 2))
12×12 SparseArrays.SparseMatrixCSC{Float64, Int64} with 52 stored entries:
  6.0  -1.0    ⋅   -1.0    ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅     ⋅
 -1.0   6.0  -1.0    ⋅   -1.0    ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅
   ⋅   -1.0   6.0    ⋅     ⋅   -1.0    ⋅     ⋅   -1.0    ⋅     ⋅     ⋅
 -1.0    ⋅     ⋅    6.0  -1.0    ⋅     ⋅     ⋅     ⋅   -1.0    ⋅     ⋅
   ⋅   -1.0    ⋅   -1.0   6.0  -1.0    ⋅     ⋅     ⋅     ⋅   -1.0    ⋅
   ⋅     ⋅   -1.0    ⋅   -1.0   6.0    ⋅     ⋅     ⋅     ⋅     ⋅   -1.0
 -1.0    ⋅     ⋅     ⋅     ⋅     ⋅    6.0  -1.0    ⋅   -1.0    ⋅     ⋅
   ⋅   -1.0    ⋅     ⋅     ⋅     ⋅   -1.0   6.0  -1.0    ⋅   -1.0    ⋅
   ⋅     ⋅   -1.0    ⋅     ⋅     ⋅     ⋅   -1.0   6.0    ⋅     ⋅   -1.0
   ⋅     ⋅     ⋅   -1.0    ⋅     ⋅   -1.0    ⋅     ⋅    6.0  -1.0    ⋅
   ⋅     ⋅     ⋅     ⋅   -1.0    ⋅     ⋅   -1.0    ⋅   -1.0   6.0  -1.0
   ⋅     ⋅     ⋅     ⋅     ⋅   -1.0    ⋅     ⋅   -1.0    ⋅   -1.0   6.0
```

"""
laplacian

laplacian(n::NTuple{N}) where {N} =
    laplacian(ntuple(i -> dirichlet, Val(N)), n)

laplacian((bc,), (n,)) = laplacian(bc, n)

function laplacian(::Dirichlet, n::Int)
    spdiagm(-1 => -ones(n-1),
            0 => 2ones(n),
            1 => -ones(n-1))
end

function laplacian(bc::NTuple{2}, n::NTuple{2})
    lap = laplacian.(bc, n)
    eye = I.(n)

    kron(eye[2], lap[1]) +
    kron(lap[2], eye[1])
end

function laplacian(bc::NTuple{3}, n::NTuple{3})
    lap = laplacian.(bc, n)
    eye = I.(n)

    kron(eye[3], eye[2], lap[1]) +
    kron(eye[3], lap[2], eye[1]) +
    kron(lap[3], eye[2], eye[1])
end

