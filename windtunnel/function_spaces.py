from dolfin import *

class FunctionSpaces(object):

    @staticmethod
    def P0(mesh):
        return FunctionSpace(mesh, 'DG', 0)
               #, vfamily='DG', vdegree=0, name='P0')

    @staticmethod
    def P1(mesh):
        return FunctionSpace(mesh, 'CG', 1)
               #, vfamily='CG', vdegree=1, name='P1')

    @staticmethod
    def P1v(mesh):
        return VectorFunctionSpace(mesh, 'CG', 1)
               #, vfamily='CG', vdegree=1, name='P1v')

    @staticmethod
    def P1DG(mesh):
        return FunctionSpace(mesh, 'DG', 1)
               #, vfamily='DG', vdegree=1, name='P1DG')

    @staticmethod
    def P1DGv(mesh):
        return VectorFunctionSpace(mesh, 'DG', 1)
               #, vfamily='DG', vdegree=1, name='P1DGv')

    @staticmethod
    def P2(mesh):
        return VectorFunctionSpace(mesh, "CG", 2)
               #, vfamily='CG', vdegree=2, name='P2')

    @staticmethod
    def P2P1(mesh):
        """ Taylor-Hood Element pair """
        return [FunctionSpaces.P2(mesh), FunctionSpaces.P1(mesh)]

    @staticmethod
    def parse_function_spaces(string, mesh):
        if string == 'P2P1':
            return FunctionSpaces.P2P1(mesh)
        else:
            raise Exception('Valid discretisation is required '
                            '- why not try "P2P1"?')


if __name__ == '__main__':

    mesh = UnitSquareMesh(10, 10)
    V, Q = FunctionSpaces.P2P1(mesh)
    print V
    print Q
    
