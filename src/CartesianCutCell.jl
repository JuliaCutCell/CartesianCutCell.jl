module CartesianCutCell

using LinearAlgebra
using SparseArrays

export dirichlet, neumann, periodic

export laplacian
export forwarddiff, backwarddiff
export forwardinterp, backwardinterp

const TupleN{T,N} = NTuple{N,T}

include("boundary.jl")
include("rectangular.jl")
include("differentiation.jl")
include("interpolation.jl")

end
