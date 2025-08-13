from .tables import *
from math import sqrt, factorial
from decimal import Decimal, getcontext

class Newton(Tables):
    __doc__="""Class containing an implementation of Newton's method for 
    integer square roots."""

    def isqrt(x):
        """For very large integers, Newton's method (or the Babylonian method) 
        can be used to calculate the integer square root and avoid overflow 
        errors associated with integer to float conversions.  This method 
        iteratively refines an initial estimate until it converges to the 
        correct integer square root.

        Although provision of this method is primarily intended for the internal
        handling of `OverflowError` exceptions for calculations involving very 
        large integers, it can nevertheless be invoked by the user provided that
        the integer argument is passed as the first positional argument rather 
        than calling the method on an instance of the class, i.e.:
        
        import PyGammaRAD
        sqrt_large_integer = PyGammaRAD.Newton.isqrt(<int>)

        Arguments:
            x: Integer object; very large intergers acceptable.

        Returns:
            An integer object corresponding to the square root of the argument 
        determined via Newton's method for handling large integers.

        Examples
            (i) To return the integer square root of 2:
            > import PyGammaRAD as pg
            > pg.Newton.isqrt(2)

            (ii) To return the integer square root of 100:
            > pg.Newton.isqrt(100)

            (iii) To return the integer square root of 2000!:
            > from math import factorial
            > pg.Newton.isqrt(factorial(2000))
        """
        if x < 0:
            raise ValueError("Square root not defined for negative numbers.")
        if not x:
            return 0

        # Initial guess based on bit length
        n = int(x)
        a, b = divmod(n.bit_length(), 2)
        x = 2**(a + b)

        y = (x + n // x) // 2
        while y < x:
            x = y
            y = (x + n // x) // 2
        return x

class ClebschGordan(Newton):
    __doc__="""MEMBER FUNCTIONS BELONGING TO THIS CLASS ARE NOT INTENDED TO 
    BE DIRECTLY INVOKED BY THE USER (*).

    Class containing methods used in the calculation of Clebsch-Gordan 
    coefficients.

    Instantiate class as:
    
        CG = ClebschGordan(j1, m1, j2, m2, j, m)
    
    to evaluate the corresponding Clebsch-Gordan coefficient:

        <j1 m1 j2 m2 | j m>

    (*): For evaluation of the Clebsch-Gordan coefficient based on the 
    `Clebsch-Gordan` class methods refer to the appropriate docstring to ensure 
    correct passage of variables to the callable:

    Method      Quantity
    ------      --------
    `cg`     :  Clebsch-Gordan coefficient
    """
    
    def __init__(self,j1,m1,j2,m2,j,m):
        self.j1, self.m1 = j1, m1
        self.j2, self.m2 = j2, m2
        self.j, self.m = j,m
        
    def delta_m(self):
        """Delta function based on magnetic substate quantum mechanical numbers,
        i.e., projections along the z-axis."""
        m1, m2, m = self.m1, self.m2, self.m
        
        delta = None
        # Test to find if "m1 + m2 = m"
        if np.fabs((m1 + m2) - m) < np.finfo(np.float32).eps: 
            delta = 1
        else: 
            delta = 0
            print("ERROR: m1 + m2 != m")
        #print("delta_m = {0}".format(delta))
        return delta

    def delta_j(self):
        """Triangular delta factor needed in the evaluation of the 
        Clebsch-Gordan coefficient."""
        j1, j2, j = self.j1, self.j2, self.j
        
        numerator = factorial((j1+j2)-j) * factorial((j1-j2)+j) * factorial((-j1)+j2+j)
        denominator = factorial(j1+j2+j+1)
        delta = sqrt(numerator/denominator)
        #print("delta_j = {0}".format(delta))
        return delta

    def coeff(self):
        """Coupling scheme function needed in the evaluation of the 
        Clebsch-Gordan coefficient."""
        j1, m1 = self.j1, self.m1
        j2, m2 = self.j2, self.m2
        j, m = self.j, self.m

        try:
        
            j1m1 = sqrt(factorial(j1+m1)*factorial(j1-m1))
            j2m2_jm = sqrt(factorial(j2+m2) * factorial(j2-m2) * factorial(j+m) * factorial(j-m) * ((2*j)+1))
            #print("coeff = {0}".format(j1m1*j2m2_jm))

        except OverflowError:
            #print("Using Newton's method")
            getcontext().prec = 1000
            j1m1 = Newton.isqrt(factorial(j1+m1)*factorial(j1-m1))
            j2m2_jm = Newton.isqrt(factorial(j2+m2) * factorial(j2-m2) * factorial(j+m) * factorial(j-m) * ((2*j)+1))
            #print("coeff = {0}".format(j1m1*j2m2_jm))
            
        return j1m1*j2m2_jm

    def run_v(self):
        """Function to return a tuple corresponding to the summation limits 
        needed for the angular momentum coupling scheme."""
        j1, j2, j = self.j1, self.j2, self.j
        m1, m2 = self.m1, self.m2
    
        start_0 = 0
        start_1 = -(j-j2+m1)
        start_2 = -(j-j1-m2)
    
        stop_0 = j1+j2-j
        stop_1 = j1-m1
        stop_2 = j2+m2
    
        start = max(start_0, max(start_1,start_2))
        stop = min(stop_0, min(stop_1,stop_2))
    
        return (start, stop)

    def sum_coupling(self):
        """Coupling scheme summation in the evaluation of the Clebsch-Gordan 
        coefficient."""
        j1, j2, j = self.j1, self.j2, self.j
        m1, m2 = self.m1, self.m2
    
        v = ClebschGordan.run_v(self)[0]
        stop = ClebschGordan.run_v(self)[1]
        #print("Run summation over integers starting at: v={0}; ending at: v={1}".format(v,int(stop)))
        sum_couple_j = 0.0
        while True:
            #print(v, sum_couple_j)
            try:
                numerator_A = (-1)**v
                denominator_A = factorial(v) * factorial(((j1+j2)-j)-v) * factorial((j1-m1)-v) * factorial((j2+m2)-v) 
        
                numerator_B = 1
                denominator_B = factorial((j-j2)+m1+v) * factorial(((j-j1)-m2)+v)
        
                sum_couple_j += (numerator_A/denominator_A) * (numerator_B/denominator_B)
                #print(v, sum_couple_j)
                v += 1
            except ValueError:
                # break loop as soon as an argument in a factorial becomes a negative value
                #print('v = {0} generates negative argument in factorial'.format(v))
                break
        
        #print("summation = {0}".format(sum_couple_j))
        return sum_couple_j
        
    def cg_calc(self):
        """Evaluate Clebsch-Gordan coefficient <j1 m1 j2 m2 | j m>"""
        j1, m1 = self.j1, self.m1
        j2, m2 = self.j2, self.m2
        j, m = self.j, self.m
        
        dm = ClebschGordan.delta_m(self)
        dj = ClebschGordan.delta_j(self)
        c = ClebschGordan.coeff(self)
        s = ClebschGordan.sum_coupling(self)

        try:
            cg = dm*dj*c*s
        except OverflowError:
            getcontext().prec = 1000
            cg = float(Decimal(dm)*Decimal(dj)*Decimal(c)*Decimal(s))
        #print("<{0} {1} {2} {3} | {4} {5}> = {6}".format("%.1f"%j1,"%.1f"%m1,"%.1f"%j2,"%.1f"%m2,"%.1f"%j,"%.1f"%m,cg))
        return cg
    
class Wigner3j(ClebschGordan):
    __doc__="""MEMBER FUNCTIONS BELONGING TO THIS CLASS ARE NOT INTENDED TO 
    BE DIRECTLY INVOKED BY THE USER (*).

    Class to handle Wigner 3j-symbols.

    Instantiate class as:
    
        W = Wigner(j1, j2, j, m1, m2, m)
        
    to evaluate the 3j symbol:
    
        (j1 j2 j
         m1 m2 m)

    (*): For evaluation of the Wigner 3-j symbol based on the `Wigner3j` class 
    methods refer to the appropriate docstring to ensure correct passage of 
    variables to the callable:

    Method      Quantity
    ------      --------
    `symb3j` :  Wigner 3-j symbol
    """
    
    def __init__(self,j1,j2,j,m1,m2,m):
        self.j1, self.j2, self.j = j1, j2, j
        self.m1, self.m2, self.m = -m1, -m2, m
        
    def symbol_3j(self):
        """Evaluate Wigner 3j-symbol."""
        j1, j2, j = self.j1, self.j2, self.j
        m1, m2, m = self.m1, self.m2, self.m
        
        cg = ClebschGordan.cg_calc(self)
        p = (-1)**(j+m+(2*j1))
        c = 1/sqrt((2*j)+1)

        try:
            symb_3j = p*c*cg
        except TypeError:
            getcontext().prec = 100000000000
            symb_3j = float(Decimal(p)*Decimal(c)*Decimal(cg))
            
        #print("({0} {1} {2}".format("%.1f"%j1, "%.1f"%j2, "%.1f"%j))
        #print(" {0} {1} {2}) = {3}".format("%.1f"%(-m1), "%.1f"%(-m2), "%.1f"%m,symb_3j))
        return symb_3j

class Racah(Wigner3j):
    __doc__ = """MEMBER FUNCTIONS BELONGING TO THIS CLASS ARE NOT INTENDED TO 
    BE DIRECTLY INVOKED BY THE USER (*).

    Class to handle Racah recoupling coeficients and Wigner 6-j symbols.
    
    Instantiate class as:
    
        W = Racah(j1, j2, j3, j4, j5, j6)
        
    to evaluate the 6j symbol:
    
        {j1 j2 j3
         j4 j5 j6}

    or the Racach coefficient:

    W(j1 j2 j5 j4; j3 j6)

    (*): For evaluation of angular momentum coefficients and symbols based on 
    the `Racah` class methods refer to the appropriate docstring to ensure 
    correct passage of variables to the callable:

    Method      Quantity
    ------      --------
    `racah`  :  Racah coefficient
    `symb6j` :  Wigner 6-j symbol
    """
    
    def __init__(self,j1,j2,j3,j4,j5,j6):
        self.j1, self.j2 = j1, j2
        self.j3, self.j4 = j3, j4
        self.j5, self.j6 = j5, j6
        
    def tri_factor(a,b,c):
        """The triangular factor (or triangular delta) Delta(j1,j2,j3) 
        provides a mathematical check to ensure given angular momentum values
        are physically compatible for forming a coupled state.  It directly 
        evaluates the validity of an angular momentum triad: Adding two angular
        momenta j1 and j2 to obtain a resultant total angular momentum j3 
        results in the formation of a (j1,j2,j3) triad.  The formation of this
        triad is subject to the triangle inequalities condition for valid 
        coupling whereupon the magnitudes must satisfy: 

        |j1 - j2| <= j3 <= j1 + j2.

        The triad thus represents the selection rules governing allowed values 
        for the total angular momentum arising from the coupling of two 
        individual angular momenta."""
        numerator = factorial((a+b)-c) * factorial((a-b)+c) * factorial((-a)+b+c)
        denominator = factorial(int(a+b+c+1))
        delta = sqrt(numerator/denominator)
        return delta
        
    def delta_product(self):
        """Evaluate product of triangular factors involved in angular momentum 
        coupling scheme."""
        delta_j1j2j3 = Racah.tri_factor(self.j1, self.j2, self.j3)
        delta_j5j4j3 = Racah.tri_factor(self.j5, self.j4, self.j3)
        delta_j1j5j6 = Racah.tri_factor(self.j1, self.j5, self.j6)
        delta_j2j4j6 = Racah.tri_factor(self.j2, self.j4, self.j6)
        delta_J = delta_j1j2j3 * delta_j5j4j3 * delta_j1j5j6 * delta_j2j4j6

        return delta_J
    
    def w(self):
        """Evaluation of quantity needed in the determination of the Racah 
        coefficient."""
        a1 = self.j1 + self.j2 + self.j3
        a2 = self.j5 + self.j4 + self.j3
        a3 = self.j1 + self.j5 + self.j6
        a4 = self.j2 + self.j4 + self.j6

        a_list = [a1,a2,a3,a4]
        max_a = max(a_list)
        
        b1 = self.j1 + self.j2 + self.j5 + self.j4
        b2 = self.j1 + self.j4 + self.j3 + self.j6
        b3 = self.j2 + self.j5 + self.j3 + self.j6

        b_list = [b1,b2,b3]
        min_b = min(b_list)
        
        z_min = int(max_a)
        z_max = int(min_b)
        
        w_coeff = 0
        OVERFLOW = False
        for z in range(z_min, z_max+1, 1):
            try:
                numerator = (-1)**(z+b1) * factorial(z+1)
                denominator = factorial(z-a1)*factorial(z-a2)*factorial(z-a3)*factorial(z-a4)*factorial(b1-z)*factorial(b2-z)*factorial(b3-z)

                ratio = numerator/denominator
                w_coeff += ratio
            except OverflowError:
                OVERFLOW = True
                getcontext().prec = 1000
                numerator = Decimal((-1)**(z+b1)) * Decimal(factorial(z+1))
                denominator = Decimal(factorial(z-a1))*Decimal(factorial(z-a2))*Decimal(factorial(z-a3))*Decimal(factorial(z-a4))*Decimal(factorial(b1-z))*Decimal(factorial(b2-z))*Decimal(factorial(b3-z))
                
                ratio = numerator/denominator
                ratio = float(ratio)
                w_coeff += ratio

        #if OVERFLOW == True:
            #print("Overflow exception handled")
                
        return w_coeff
    
    def phase(self):
        """Evaluate phase factor."""
        return (-1)**(self.j1+self.j2+self.j4+self.j5)
    
    def W(self):
        """Evaluate Racah coefficient."""
        return Racah.delta_product(self)*Racah.w(self)
    
    def symbol_6j(self):
        """Evaluate Wigner 6-j symbol."""
        return Racah.W(self)*Racah.phase(self)

class Wigner9j(Racah):
    __doc__ = """MEMBER FUNCTIONS BELONGING TO THIS CLASS ARE NOT INTENDED TO 
    BE DIRECTLY INVOKED BY THE USER (*).

    Class to handle Wigner 9-j symbol.

    Instantiate class as:
    
        W = Wigner9j(j1, j2, j3, j4, j5, j6, j7, j8, j9)
        
    to evaluate the 9j symbol:
    
        {j1 j2 j3
         j4 j5 j6
         j7 j8 j9}

    (*): For evaluation of the Wigner 9-j symbol based on the `Wigner9j` class 
    methods refer to the appropriate docstring to ensure correct passage of 
    variables to the callable:

    Method      Quantity
    ------      --------
    `symb9j` :  Wigner 9-j symbol
    """
    
    def __init__(self,j1,j2,j3,j4,j5,j6,j7,j8,j9):
        self.j1, self.j2, self.j3 = j1, j2, j3
        self.j4, self.j5, self.j6 = j4, j5, j6 
        self.j7, self.j8, self.j9 = j7, j8, j9

    def symbol_9j(self):
        """Evaluate Wigner 9-j symbol."""
        imax = int(min(self.j1+self.j9, self.j2+self.j6, self.j4+self.j8) * 2)
        imin = imax % 2
        sum_res = 0
        #print(imin, imax)
        for x in range(int(imin), int(imax)+1, 2):
            try:
                W1 = Racah(self.j1, self.j4, self.j7, self.j8, self.j9, x/2)
                W2 = Racah(self.j2, self.j5, self.j8, self.j4, x/2, self.j6)
                W3 = Racah(self.j3, self.j6, self.j9, x/2, self.j1, self.j2)
            
                sum_res = sum_res + (x+1) * W1.symbol_6j() * W2.symbol_6j() * W3.symbol_6j()
                
            except ValueError:
                pass
            
        return sum_res
