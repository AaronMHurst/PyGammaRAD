from .tables import *
from .am_formulae import *

class AngularMomentumCalculations(Tables):
    __doc__="""Class containing methods for coupling and recoupling of angular 
    momenta."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)        

    # Implement basic AM calculations here:
    def cg(self,j1,m1,j2,m2,j,m):
        """The Clebsch-Gordan coefficient used to combine different angular 
        momenta may be evaluated by entering all j and m terms in the 
        order in which they appear in the corresponding coefficient.

        To evaluate: <j1 m1 j2 m2 | j3 m3> 

        Call the method as: `cg(j1,m1,j2,m2,j3,m3)`

        The above Clebsch-Gordan coefficient is closely related to the 
        Wigner 3-j symbol:

        (j1 j2 j3
         -m1 -m2 m3)

        Arguments:
            ji: Set of three angular momentum quantum numbers used for coupling.
                Numerical data types should be entered as floats or integers.
            mi: Set of three corresponding angular momentum projections used 
                for coupling (magnetic quantum numbers). Numerical data types 
                should be entered as floats or integers.

        Returns:
            The Clebsch-Gordan coefficient as a floating-point object.

        Example:
            To evaluate the  Clebsch-Gordan coefficient 
            <j1=5/2 m1=3/2 j2=5/2 m2=-1/2 | j3=1 m3=1>:
        
            > cg(2.5, 1.5, 2.5, -0.5, 1, 1)
        """
        self.j1, self.m1 = j1, m1
        self.j2, self.m2 = j2, m2
        self.j, self.m = j, m

        CG = ClebschGordan(self.j1, self.m1, self.j2, self.m2, self.j, self.m)
        return CG.cg_calc()

    def symb3j(self,j1,j2,j,m1,m2,m):
        """The Wigner 3-j symbol used to calculate the coupling of different 
        angular momenta may be evaluated by entering all j and m terms in the 
        order in which they appear in the corresponding 3-j symbol.

        To evaluate: (j1 j2 j3 
                      m1 m2 m3)

        Call the method as: `symb3j(j1,j2,j3,m1,m2,m3)`

        The above 3-j symbol is closely related to the Clebsch-Gordan 
        coefficient:

        <j1 -m1 j2 -m2|j3 m3>

        Arguments:
            ji: Set of three angular momentum quantum numbers used for coupling.
                Numerical data types should be entered as floats or integers.
            mi: Set of three corresponding angular momentum projections used 
                for coupling (magnetic quantum numbers). Numerical data types 
                should be entered as floats or integers.

        Returns:
            The Wigner 3-j symbol as a floating-point object.

        Example:
            To evaluate the 3-j symbol {j1=7/2 j2=5/2 j3=2 
                                        m1=7/2 m2=-3/2 m3=-2}:
        
            > symb3j(3.5, 2.5, 2, 3.5, -1.5, -2)

        """
        self.j1, self.j2, self.j = j1, j2, j
        self.m1, self.m2, self.m = m1, m2, m

        W = Wigner3j(self.j1, self.j2, self.j, self.m1, self.m2, self.m)
        return W.symbol_3j()

    def racah(self,j1,j2,j3,j4,j5,j6):
        """The Racah recoupling coefficient describes the transformation 
        between different angular momenta coupling schemes and may be 
        evaluated by entering all terms in the order in which they appear in 
        Racah-W coefficient.

        To evaluate: W(j1 j2 j3 j4; j5 j6)

        Call the method as: `racach(j1,j2,j3,j4,j5,j6)`

        The above Racah coefficient is closely related to the 6-j symbol:
        {j1 j2 j5
         j4 j3 j6}.

        Arguments:
            ji: Set of six angular momenta used for coupling.  Numerical data 
                types should be entered as floats or integers.

        Returns:
            The Racach recoupling coefficient as a floating-point object.

        Example:
            To evaluate the Racah W(j1=15 j2=15 j3=17 j4=17; j5=2 j6=3):
        
            > racah(15, 15, 17, 17, 2, 3)
        """
        self.j1, self.j2, self.j3 = j1, j2, j3
        self.j4, self.j5, self.j6 = j4, j5, j6

        W = Racah(self.j1, self.j2, self.j5, self.j4, self.j3, self.j6)
        return W.W()

    def symb6j(self,j1,j2,j3,j4,j5,j6):
        """The Wigner 6-j symbol based on the combination of different angular 
        momenta coupling schemes may be evaluated by entering all terms in the 
        order in which they appear in the corresponding 6-j symbol.

        To evaluate: {j1 j2 j3 
                      j4 j5 j6}

        Call the method as: `symb6j(j1,j2,j3,j4,j5,j6)`

        The above 6-j symbol is closely related to the Racah coefficient:
        W(j1 j2 j5 j4; j3 j6).

        Arguments:
            ji: Set of six angular momenta used for coupling.  Numerical data 
                types should be entered as floats or integers.

        Returns:
            The Wigner 6-j symbol as a floating-point object.

        Example:
            To evaluate the 6-j symbol {j1=6 j2=6 j3=4 
                                        j4=9/2 j5=7/2 j6=11/2}:
        
            > symb6j(6, 6, 4, 4.5, 3.5, 5.5)

        """
        self.j1, self.j2, self.j3 = j1, j2, j3
        self.j4, self.j5, self.j6 = j4, j5, j6

        W = Racah(self.j1, self.j2, self.j3, self.j4, self.j5, self.j6)
        return W.symbol_6j()

    def symb9j(self,j1,j2,j3,j4,j5,j6,j7,j8,j9):
        """The Wigner 9-j symbol based on the combination of different angular 
        momenta coupling schemes may be evaluated by entering all terms in the 
        order in which they appear in the corresponding 9-j symbol.

        To evaluate: {j1 j2 j3 
                      j4 j5 j6
                      j7 j8 j9}

        Call the method as: `symb9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)`

        Arguments:
            ji: Set of nine angular momenta used for coupling.  Numerical data 
                types should be entered as floats or integers.

        Returns:
            The Wigner 9-j symbol as a floating-point object.

        Example:
            To evaluate the 9-j symbol {j1=3 j2=4 j3=2 
                                        j4=7/2 j5=7/2 j6=2
                                        j7=1/2 j8=1/2 j9=1}:
        
            > symb6j(3, 4, 2, 3.5, 3.5, 2, 0.5, 0.5, 1)
        """
        self.j1, self.j2, self.j3 = j1, j2, j3
        self.j4, self.j5, self.j6 = j4, j5, j6
        self.j7, self.j8, self.j9 = j7, j8, j9

        W = Wigner9j(self.j1, self.j2, self.j3, self.j4, self.j5, self.j6, self.j7, self.j8, self.j9)
        return W.symbol_9j()

    
