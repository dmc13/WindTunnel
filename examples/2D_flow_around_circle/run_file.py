from windtunnel import *

# Grab default parameter set
problem_parameters = Problem.default_parameters()

# Set up the domain
domain = FileDomain('mesh/mesh.xml')
problem_parameters.domain = domain

# Choose our discretisation 
problem_parameters.discretisation = "P2P1"

# Set boundary conditions, vel; no-slip on boundarys, press; inflow then outflow
bcs = BoundaryConditions()
bcs.add_bc_u((0, 0), facet_id=3)
p_in = Expression("sin(3.0*t)", t=0.0)
bcs.add_bc_p(p_in, facet_id=1, time_dependent=True)
bcs.add_bc_p(0, facet_id=2)
problem_parameters.bcs = bcs

# Set up other parameters
problem_parameters.dt = 0.01
problem_parameters.finish_time = 3 
problem_parameters.viscosity = 0.01

# Set up problem
problem = Problem(problem_parameters)

# Grab default solver parameter set
solver_parameters = Solver.default_parameters()
solver_parameters.live_plotting = True
solver_parameters.dump_period = 1000

# Set up solver
solver = Solver(solver_parameters, problem)

# Solve
solver.solve()

