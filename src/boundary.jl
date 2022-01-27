abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition end
struct Periodic <: BoundaryCondition end

const dirichlet = Dirichlet()
const periodic = Periodic()

has(::Dirichlet) = true
has(::Periodic) = false

