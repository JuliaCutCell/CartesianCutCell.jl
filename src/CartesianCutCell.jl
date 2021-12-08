module CartesianCutCell

const TupleN{T,N} = NTuple{N,T}

export laplacian

include("boundary.jl")
include("rectangular.jl")

end
