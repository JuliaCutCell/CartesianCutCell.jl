struct Stencil{W} end

width(::Stencil{W}) where {W} = W

# Marker and Cell
const mac = Stencil{1}()

stencil(::NTuple{N,Any}) where {N} = ntuple(i -> mac, Val(N))

