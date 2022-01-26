module CartesianCutCell

using LinearAlgebra
using SparseArrays
using FillArrays
using LazyArrays

export dirichlet, periodic

export mesh
#export mask
export pad

export laplacian
export forwarddiff, backwarddiff
export forwardinterp, backwardinterp

include("boundary.jl")
include("mesh.jl")
#include("mask.jl")
include("pad.jl")
include("rectangular.jl")
include("differentiation.jl")
include("interpolation.jl")

end
