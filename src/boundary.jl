abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition
    width::Int
end
struct Periodic <: BoundaryCondition end

const dirichlet = Dirichlet(1)
const periodic = Periodic()

Base.length(bc::Dirichlet) = getproperty(bc, :width)

