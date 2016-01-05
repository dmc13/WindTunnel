from dolfin import *
from helpers import FrozenClass

class ProblemParameters(FrozenClass):
    """
    """
    domain = None
    finish_time = 3
    dt = 0.01
    viscosity = 0.01
    bcs = None

class Problem(object):
    """
    """

    def __init__(self, parameters):
	
	if not isinstance(parameters, ProblemParameters):
	    raise TypeError, "Problem requires parameters of \
			     type ProblemParameters"
        else:
            self.parameters = parameters
        
        # Generate the boundary conditions
        self.parameters.bcs.generate_bcs()
    
    @staticmethod
    def default_parameters():
	return ProblemParameters()


