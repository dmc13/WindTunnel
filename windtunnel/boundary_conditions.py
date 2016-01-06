from dolfin import *
from function_spaces import FunctionSpaces

class BoundaryConditions(object):

    def __init__(self):
        self.bc_u_list = []
        self.bc_p_list = []
        self.domain = None 
        self.function_spaces = None

    def add_bc_u(self, expression=None, facet_id=None, time_dependent=False):
        # Velocity bcs
        #bc = DirichletBC(self.V, expression, facet_id)
        self.bc_u_list.append([expression, facet_id, time_dependent])

    def add_bc_p(self, expression=None, facet_id=None, time_dependent=False):
        # Pressure bcs
        #bc = DirichletBC(self.Q, expression, facet_id)
        self.bc_p_list.append([expression, facet_id, time_dependent])

    def generate_bcs(self):
        self.bc_u, self.bc_p = [], []
        V, Q = self.function_spaces
        for i in range(len(self.bc_u_list)):
            self.bc_u.append(DirichletBC(V, self.bc_u_list[i][0],
                                         self.domain.facet_ids, self.bc_u_list[i][1]))
        for i in range(len(self.bc_p_list)):
            self.bc_p.append(DirichletBC(Q, self.bc_p_list[i][0],
                                         self.domain.facet_ids, self.bc_p_list[i][1]))

    def update_time(self, t):
        # Update the time in the expression for temporal BCs
        # TODO this seems like a poor way of doing this...
        for bc in self.bc_u_list:
            if bc[2]:
                bc[0].t = t
        for bc in self.bc_p_list:
            if bc[2]:
                bc[0].t = t
        self.generate_bcs
#        from IPython import embed; embed()

