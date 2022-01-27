function integrate(func, x, dx)
    V = map(*, dx...)
    A = _integrate(func, x, dx)
    V, A 
end

_integrate(func, x, dx::NTuple{1,Any}) = (Ones{Bool}(dx[1]),)

_integrate(func, x, dx::NTuple{2,Any}) = (dx[2], dx[1])

_integrate(func, x, dx::NTuple{3,Any}) =
    (dx[2] .* dx[3], dx[3] .* dx[1], dx[1] .* dx[2])

