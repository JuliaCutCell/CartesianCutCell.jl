module CartesianCutCell

using LinearAlgebra
using SparseArrays
using FillArrays
using LazyArrays

export dirichlet, periodic
export mac, stencil

export mesh
export restrict

export laplacian
export forwarddiff, backwarddiff
export forwardinterp, backwardinterp

include("boundary.jl")
include("stencil.jl")
include("mesh.jl")
include("differentiation.jl")
include("interpolation.jl")
include("restrict.jl")

end
