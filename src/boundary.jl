abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition end
struct Neumann <: BoundaryCondition end
struct Periodic <: BoundaryCondition end

const dirichlet = Dirichlet()
const neumann = Neumann()
const periodic = Periodic()

