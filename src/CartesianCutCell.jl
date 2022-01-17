module CartesianCutCell

using LinearAlgebra
using SparseArrays
using LazyArrays

export dirichlet, neumann, periodic

export mesh

export laplacian
export forwarddiff, backwarddiff
export forwardinterp, backwardinterp

const TupleN{T,N} = NTuple{N,T}

include("boundary.jl")
include("grid.jl")
include("rectangular.jl")
include("differentiation.jl")
include("interpolation.jl")

end
