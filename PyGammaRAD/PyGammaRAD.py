from .tables import *
from .am_formulae import *
from .am_methods import *
from .angular_distributions import *
from .log_handlers import *

class AngularMomentum(Legendre):    
    __doc__="""Class to handle the coupling and recoupling schemes of angular 
    momenta needed in the calculation of gamma-ray angular distributions in 
    aligned nuclei.  The angular momentum calculators also have purpose in 
    wider quantum mechanical applications."""

    def __init__(self):
        Yamazaki.__init__(self)
        RoseAndBrink.__init__(self)
        AngularMomentumCalculations.__init__(self)
        AngularDistributions.__init__(self)
        Legendre.__init__(self)

