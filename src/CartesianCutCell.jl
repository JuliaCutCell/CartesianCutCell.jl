module CartesianCutCell

using SparseArrays

export laplacian

const TupleN{T,N} = NTuple{N,T}

include("boundary.jl")
include("rectangular.jl")

end
