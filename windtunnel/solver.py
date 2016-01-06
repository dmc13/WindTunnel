from dolfin import *
import os
from helpers import FrozenClass, StateWriter
from problem import Problem
from function_spaces import FunctionSpaces

class SolverParameters(FrozenClass):
    
    # Some fenics housekeeping:
    # - keep the std_out tidy when in parallel
    parameters["std_out_all_processes"] = False
    # - use amg preconditioner if available
    prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

    dump_period = 1
    output_dir = os.curdir 
    live_plotting = False


class Solver(object):

    def __init__(self, parameters, problem):

        # Check inputs and attach to class instance
	if not isinstance(parameters, SolverParameters):
	    raise TypeError, "Solver requires parameters of \
			     type SolverParameters"
        else: 
	    self.parameters = parameters
	if not isinstance(problem, Problem):
	    raise TypeError, "Problem must be of type Problem"
        else:
            self.problem = problem

        # Make some things more accessible:
        self.bcs = self.problem.parameters.bcs


    @staticmethod
    def default_parameters():
        return SolverParameters()

    
    def solve(self):

        # Fetch the function spaces
        V, Q = self.problem.parameters.discretisation
        # Set up the velocity functions 
        u = TrialFunction(V)
        v = TestFunction(V)
        u0 = Function(V)
        u1 = Function(V) 
	# Set up the pressure functions
        p = TrialFunction(Q)
        q = TestFunction(Q)
        p1 = Function(Q)

        # Define coefficients
	k = Constant(self.problem.parameters.dt)
	f = Constant((0,0))
	nu = Constant(self.problem.parameters.viscosity)

	# Tentative velocity step
	F1 = (1/k) * inner(u-u0, v)*dx          \
	     + inner(grad(u0)*u0, v)*dx         \
	     + nu * inner(grad(u), grad(v))*dx  \
	     - inner(f, v)*dx
	a1 = lhs(F1)
	L1 = rhs(F1)
	A1 = assemble(a1)

	# Pressure update
	a2 = inner(grad(p), grad(q))*dx
	L2 = (-1./k) * div(u1) * q * dx
	A2 = assemble(a2)

	# Velocity update
	a3 = inner(u, v)*dx
	L3 = inner(u1, v)*dx - k * inner(grad(p1), v)*dx
	A3 = assemble(a3)

        # Iterate over timesteps, start at dt
        self.t = self.problem.parameters.dt
        while self.t < self.problem.parameters.finish_time:
            print 'Starting timestep at time = ', self.t

            # Update the boundary conditions for the current timestep
            self.bcs.update_time(t=self.t)

            # Tentative velocity step
            print 'Step 1: Computing tentative velocity'
            b1 = assemble(L1)
            [bc.apply(A1, b1) for bc in self.bcs.bc_u]
            solve(A1, u1.vector(), b1, "gmres", "default")
            end()

            # Pressure update
            print 'Step 2: Updating pressure'
            b2 = assemble(L2)
            [bc.apply(A2, b2) for bc in self.bcs.bc_p]
            solve(A2, p1.vector(), b2, "gmres", self.parameters.prec)
            end()

	    # Velocity update
            print 'Step 3: Updating velocity'
            b3 = assemble(L3)
            [bc.apply(A3, b3) for bc in self.bcs.bc_u]
            solve(A3, u1.vector(), b3, "gmres", "default")
            end()

            # Timestep update
            print 'Step 4: Dump and update timestep'
            if (self.t/self.problem.parameters.dt)%self.parameters.dump_period == 0:
                writer = StateWriter(solver=self)
                writer.write(u1, p1)
            u0.assign(u1)
            self.t += self.problem.parameters.dt
            
            if self.parameters.live_plotting:
                plot(p1, title="Pressure", rescale=True)
                plot(u1, title="Velocity", rescale=True)

        if self.parameters.live_plotting:
            interactive()



