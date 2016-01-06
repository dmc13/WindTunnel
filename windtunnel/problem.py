from dolfin import *
from helpers import FrozenClass
from function_spaces import FunctionSpaces

class ProblemParameters(FrozenClass):
    """
    """
    domain = None
    discretisation = None
    finish_time = 3
    dt = 0.01
    viscosity = 0.01
    bcs = None

class Problem(object):
    """
    """

    def __init__(self, parameters):
	
	if not isinstance(parameters, ProblemParameters):
	    raise TypeError, "Problem requires parameters of " \
			     "type ProblemParameters"
        else:
            self.parameters = parameters

        if self.parameters.domain == None:
            raise NotImplementedError, "Problem parameters requires a domain to "\
                                       "be defined"
        if self.parameters.discretisation == None:
            raise NotImplementedError, "Problem parameters requires "\
                                       "function_spaces to be defined"
        if self.parameters.bcs == None:
            raise NotImplementedError, "Problem parameters requires a bcs to be "\
                                       "defined"
        
        # Parse the function spaces string input
        self.parameters.discretisation = \
            FunctionSpaces.parse_function_spaces(self.parameters.discretisation,
                                                 self.parameters.domain.mesh)
        
        # Generate the boundary conditions
        self.parameters.bcs.domain = self.parameters.domain
        self.parameters.bcs.function_spaces = self.parameters.discretisation
        self.parameters.bcs.generate_bcs()

    
    @staticmethod
    def default_parameters():
	return ProblemParameters()


