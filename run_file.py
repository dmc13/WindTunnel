from daves_navier_stokes import *

# keep the std_out tidy when in parallel
parameters["std_out_all_processes"] = False

# mesh
mesh = Mesh("lshape.xml.gz")


