from .PyGammaRAD import *
from .tables import *
from .am_formulae import *
from .am_methods import *
from .angular_distributions import *

__version__='0.1.0'
__author__='Aaron M. Hurst'

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """Method to return absolute path of the data files inside the root of 
    the Python package."""
    return os.path.join(_ROOT, path)
