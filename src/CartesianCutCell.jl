module CartesianCutCell

using LinearAlgebra
using SparseArrays
using FillArrays
using LazyArrays

export dirichlet, neumann, periodic

export mesh
export mask
export pad

export laplacian
export forwarddiff, backwarddiff
export forwardinterp, backwardinterp

const TupleN{T,N} = NTuple{N,T}

include("boundary.jl")
include("grid.jl")
include("mask.jl")
include("pad.jl")
include("rectangular.jl")
include("differentiation.jl")
include("interpolation.jl")

end
