from .PyGammaRAD import *
from .tables import *
from .am_formulae import *
from .am_methods import *
from .angular_distributions import *
from .log_handlers import *

__version__='0.2.0'
__author__='Aaron M. Hurst'

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """Method to return absolute path of the data files inside the root of 
    the Python package."""
    return os.path.join(_ROOT, path)

logger.info(f"\n-------------------------------------------------------\n Welcome to PyGammaRAD v.{__version__} \n\n To write logs to file during session call the method:\n\n PyGammaRAD.write_logs() \n\n For an overview of the project visit the repo: \n\n https://github.com/AaronMHurst/PyGammaRAD \n-------------------------------------------------------")
