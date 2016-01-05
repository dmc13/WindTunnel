from windtunnel import *

# Grab default parameter set
problem_parameters = Problem.default_parameters()

# Set up the domain
domain = FileDomain('mesh/mesh.xml')
problem_parameters.domain = domain

# Choose our discretisation
function_spaces = FunctionSpaces.P2P1(domain.mesh)

# Set boundary conditions, vel; no-slip on boundaries, press; inflow then outflow
bcs = BoundaryConditions(function_spaces, domain)

bcs.add_bc_u((0, 0), facet_id=3)

p_in = Expression("((-6*x[1]*(-0.41*x[1]))*sin((pi*t)/8.0)/0.02825761)", t=0.0)

bcs.add_bc_p(p_in, facet_id=1, time_dependent=True)
bcs.add_bc_p(0, facet_id=2)

problem_parameters.bcs = bcs

# Set up other parameters
problem_parameters.dt = 0.01
problem_parameters.finish_time = 8 
problem_parameters.viscosity = 1./1000.

# Set up problem
problem = Problem(problem_parameters)

# Grab default solver parameter set
solver_parameters = Solver.default_parameters()

# Set up solver
solver = Solver(solver_parameters, problem)

# Solve
solver.solve()

