using LinearAlgebra
using SparseArrays

function laplacian()
    bc = (dirichlet, periodic)
    st = stencil(bc)
    n = (3, 4)

    n̅, x̅, x̅ₘ, dx̅ = mesh(bc, st, n)

    δ̅₋, δ̅₊ = backwarddiff(bc, n̅), forwarddiff(bc, n̅)
    σ̅₋, σ̅₊ = backwardinterp(bc, n̅), forwardinterp(bc, n̅)

    Δ̅ = -mapreduce(*, +, δ̅₋, δ̅₊)
    Δ = restrict(Δ̅, n̅, n) |> sparse

    x = restrict.(x̅, Ref(n̅), Ref(n))
    xₘ = restrict.(x̅ₘ, Ref(n̅), Ref(n))

    f(u, v) = π * u ^ 2 * v - √2 * v ^ 2 - √3u + 1

    y = f.(x...) |> collect
    b = Δ * y
    yₕ = Δ \ b

    norm(y - yₕ)
end

@testset "rectangular" begin
    @test isapprox(laplacian(), 0, atol=1e-15)
end

