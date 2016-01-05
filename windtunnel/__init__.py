"""

WindTunnel - a simple Navier-Stokes solver

"""

__version__ = 'alpha'
__author__  = ''
__credits__ = ['']
__license__ = 'GPL-3'
__maintainer__ = ''
__email__ = ''

from function_spaces import *
from helpers import *
from problem import *
from solver import *
from boundary_conditions import *

import numpy as np
np.random.seed(21)
from IPython import embed
from matplotlib import pyplot as plt

