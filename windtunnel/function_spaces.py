from dolfin import *

class FunctionSpaces(object):

    def __init__(self):
        # Creates the function spaces
        self.P0 = FunctionSpace(self.mesh, 'DG', 0, vfamily='DG', vdegree=0,
                                name='P0')
        self.P1 = FunctionSpace(self.mesh, 'CG', 1, vfamily='CG', vdegree=1,
                                name='P1')
        self.P1v = VectorFunctionSpace(self.mesh, 'CG', 1, vfamily='CG', vdegree=1,
                                       name='P1v')
        self.P1DG = FunctionSpace(self.mesh, 'DG', 1, vfamily='DG', vdegree=1,
                                  name='P1DG')
        self.P1DGv = VectorFunctionSpace(self.mesh, 'DG', 1, vfamily='DG', vdegree=1,
                                         name='P1DGv')
        self.P2 = VectorFunctionSpace(mesh, "CG", 2, vfamily='CG', vdegree=2,
                                      name='P2')


        self.append((function_name, expression, facet_id, bctype))
