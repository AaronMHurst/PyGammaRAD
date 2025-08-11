from .tables import *
from .am_formulae import *
from .am_methods import *
from .angular_distributions import *

#class AngularMomentum(AngularMomentumCalculations):
#class AngularMomentum(AngularDistributions):
class AngularMomentum(Legendre):    
    __doc__="""Class to handle the coupling of angular momentum vectors used in 
    the calculation of gamma-ray angular distributions in aligned nuclei."""

    def __init__(self):
        Tables.__init__(self)
        AngularMomentumCalculations.__init__(self)
        AngularDistributions.__init__(self)
        Legendre.__init__(self)
        #Newton.__init__(self)
        #ClebschGordan.__init__(self)
        #Wigner3j.__init__(self)
        #Racah.__init__(self)
        #Wigner9j.__init__(self)
        

