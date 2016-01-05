from dolfin import *
from function_spaces import FunctionSpaces

class BoundaryConditionSet(list):
    # Taken from OTF
    def add_bc(self, function_name, 
               expression=None, facet_id=None, bctype="strong_dirichlet"):
        """ Valid choices for bctype: "weak_dirichlet", "strong_dirichlet",
            "flather", "free_slip"
        """
        if expression is None and bctype!="free_slip":
            raise TypeError("Boundary condition of type %s requires "
                            "expression argument." % bctype)
        if expression is not None and bctype=="free_slip":
            raise TypeError('Boundary condition of type "free_slip" '
                            'does not allow expression argument.')
        if facet_id is None:
            raise TypeError('facet_id argument to add_bc() method is '
                            'not optional')

	self.append((function_name, expression, facet_id, bctype))

class BoundaryConditions(object):

    def __init__(self):
        function_spaces = FunctionSpaces()
        self.V = function_spaces.P2
        self.bc_u = []
        self.Q = function_spaces.P1
        self.bc_p = []

    def add_bc_u(self, expression=None, facet_id=None, bctype="strong_dirichlet"):
        # Velocity bcs
        bc = DirichletBC(self.V, expression, facet_id, bctype)
        self.bc_u.append(bc)

    def add_bc_p(self, expression=None, facet_id=None, bctype="strong_dirichlet"):
        # Pressure bcs
        bc = DirichletBC(self.Q, expression, facet_id, bctype)
        self.bc_p.append(bc)

    def update_time(self, t):
        # Update the time in the expression for temporal BCs
        # TODO not sure we can access the expressions like this...
        for bc in self.bc_u:
            if hasattr(bc[1], "t"):
                bc[1].t = t
        for bc in self.bc_p:
            if hasattr(bc[1], "t"):
                bc[1].t = t

