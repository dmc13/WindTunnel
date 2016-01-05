import os.path
import dolfin
from dolfin import * #TODO tidy up this import

class Domain(object):
    """ An abstract domain class. """

    def __init__(self):
        raise NotImplementedError("Domain is a base class only.")

    def __str__(self):
        comm = mpi_comm_world()
        hmin = MPI.min(comm, self.mesh.hmin())
        hmax = MPI.max(comm, self.mesh.hmax())
        num_cells = MPI.sum(comm, self.mesh.num_cells())

        s = "Number of mesh elements: %i.\n" % num_cells
        s += "Mesh element size: %f - %f." % (hmin, hmax)

        return s

    @property
    def ds(self):
        """A :class:`dolfin.Measure` for the facet parts of the domain."""
        return self._ds


    @property
    def dx(self):
        """A :class:`dolfin.Measure` for the cell subdomains."""
        return self._dx


class FileDomain(Domain):
    """ Create a domain from DOLFIN mesh files (.xml).
    :param mesh_file: The .xml file of the mesh.
    :type mesh_file: str
    :param facet_ids_file: The .xml file containing the facet ids of the mesh.
        If None, the default is to `mesh_file` + "_facet_region.xml".
    :type facet_ids_file: str
    :param cell_ids_file: The .xml file containing the cell ids of the mesh.
        If None, the default is to `mesh_file` + "_physical_region.xml".
    :type cell_ids_file: str
    """

    def __init__(self, mesh_file, facet_ids_file=None, cell_ids_file=None):

        #: A :class:`dolfin.Mesh` containing the mesh.
        self.mesh = dolfin.Mesh(mesh_file)

        # Read facet markers
        if facet_ids_file is None:
            facet_ids_file = (os.path.splitext(mesh_file)[0] +
                              "_facet_region.xml")

        # Read cell markers
        if cell_ids_file is None:
            cell_ids_file = (os.path.splitext(mesh_file)[0] +
                            "_physical_region.xml")

        #: A :class:`dolfin.FacetFunction` containing the surface markers.
        self.facet_ids = dolfin.MeshFunction("size_t", self.mesh, facet_ids_file)
        #: A :class:`dolfin.Measure` for the facet parts.
        self._ds = dolfin.Measure('ds')[self.facet_ids]

        #: A :class:`dolfin.CellFunction` containing the area markers.
        self.cell_ids = dolfin.MeshFunction("size_t", self.mesh, cell_ids_file)
        #: A :class:`dolfin.Measure` for the cell subdomains.
        self._dx = dolfin.Measure("dx")[self.cell_ids]


class RectangularDomain(Domain):
    """ Create a rectangular domain.
    :param x0: The x coordinate of the bottom-left.
    :type x0: float
    :param y0: The y coordinate of the bottom-left.
    :type y0: float
    :param x1: The x coordinate of the top-right corner.
    :type x1: float
    :param y1: The y coordinate of the top-right corner.
    :type y1: float
    """


    def __init__(self, x0, y0, x1, y1, nx, ny):
        #: A :class:`dolfin.Mesh` containing the mesh.
        mpi_comm = dolfin.mpi_comm_world()
        self.mesh = dolfin.RectangleMesh(mpi_comm, dolfin.Point(x0, y0),
                dolfin.Point(x1, y1), nx, ny)

        class Left(dolfin.SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and dolfin.near(x[0], x0)

        class Right(dolfin.SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and dolfin.near(x[0], x1)

        class Sides(dolfin.SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and (dolfin.near(x[1], y0) or dolfin.near(x[1], y1))

        # Initialize sub-domain instances
        left = Left()
        right = Right()
        sides = Sides()

        # Create facet markers
        #: A :class:`dolfin.FacetFunction` containing the surface markers.
        self.facet_ids = dolfin.FacetFunction('size_t', self.mesh)
        self.facet_ids.set_all(0)
        left.mark(self.facet_ids, 1)
        right.mark(self.facet_ids, 2)
        sides.mark(self.facet_ids, 3)
        #: A :class:`dolfin.Measure` for the facet parts.
        self._ds = dolfin.Measure('ds')[self.facet_ids]

        #: A :class:`dolfin.CellFunction` containing the area markers.
        self.cell_ids = dolfin.CellFunction("size_t", self.mesh)
        self.cell_ids.set_all(0)
        #: A :class:`dolfin.Measure` for the cell cell_ids.
        self._dx = dolfin.Measure("dx")[self.cell_ids]
