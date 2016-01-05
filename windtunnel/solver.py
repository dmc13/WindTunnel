from dolfin import *
from helpers import FrozenClass
from problem import Problem
from function_spaces import FunctionSpaces

class SolverParameters(FrozenClass):
    
    pass


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


    @staticmethod
    def default_parameters():
        return SolverParameters()

    
    def solve(self, parameters, problem):

        # Set up the velocity functions and spaces
        function_space = FunctionSpaces()
        V = function_space.P2
        u = TrialFunction(V)
        v = TestFunction(V)
        u0 = Function(V)
        u1 = Function(V) 
	# Set up the pressure functions and spaces
        Q = function_space.P1
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


	     



