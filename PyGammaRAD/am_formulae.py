from .tables import *
from math import sqrt, factorial, pi, lgamma, exp, log
from decimal import Decimal, getcontext
import sys
sys.set_int_max_str_digits(0) # disable limit; only limitation available memory

class Factorial(RoseAndBrink):
    __doc__="""Class containing different implementations of methods for 
    calculating the factorial of a non-negative integer.  By default, the 
    `PyGammaRAD` library uses a function that returns the native 
    `math.factorial` method.  However, this class is intended to allow users 
    to test or adjust different methods for calculating the factorial according 
    to preference.
    """
    
    def factorial_n(number):
        """That native `math.factorial` function n! is defined for 
        non-negative integers and calculates the product of all positive 
        integers less than or equal to n.  This method gets returned whenever 
        a `PyGammaRAD` method calls a factorial function and ensures the 
        argument gets intepreted as an integer."""
        return factorial(int(number))

    def factorial_recursive(number,MAX_DEPTH=2000):
        """Recursive method to calculate the factorial of non-negative integers.

        WARNING: This function may be very inefficient compared to the native 
        `math.factorial` method if the default recursion limit is exceeded.  
        Setting the maximum recursion depth too high will likely result in a 
        very long execution time.  For small recursion depth values this method 
        is probably fine, for larger values it is probably best to use the 
        native factorial method.

        Args:
            number: Number object (integer or float) for which to calculate the 
                    factorial.
            MAX_DEPTH: The recursion limit for calculating the factorial.  By 
                       default the integer value is set to 2000.

        Returns:
            The factorial of the `number` argument passed.

        Raises:
            Negative `number` arguments raise a ValueError exception.
        """
        number = int(number)
        MAX_DEPTH = int(MAX_DEPTH)
        if number < 0:
            raise ValueError("Factorial for negative numbers not defined.")
        elif number == 0 or number == 1:
            return 1
        else:
            try:
                return number * Factorial.factorial_recursive(number - 1)
            except RecursionError:
                import sys
                sys.setrecursionlimit(MAX_DEPTH)
                return number * Factorial.factorial_recursive(number - 1)

    def factorial_gamma(number):
        """The Gamma function extends the concept of a generalized factorial 
        that can handle complex numbers and non-negative floating-point 
        numbers.  Both the Gamma and factorial functions satisfy the relation: 

            Gamma(n+1) = n!
        
        This method can be used as an alternative to the factorial methods 
        when a float object gets passed that cannot be interpreted as an 
        integer.

        Args:
            number: Number object (integer or float) for which to calculate the 
                    (n+1) factorial using the `math.gamma` function.

        Returns:
            The factorial of the (`number` + 1) argument passed returned as a 
            float.

        Raises:
            Large `number` arguments (n>170) raise an OverflowError exception.
        """
        from math import gamma
        if number <= 0:
            raise ValueError("Gamma function not defined for zero or negative numbers.")
        else:
            return gamma(number + 1)

class Newton(Factorial):
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

    def data_check(x1,x2,x3,x4,x5,x6,*args):
        """Base level check of data to ensure all entries are given as integral 
        or half-integral values only, i.e., 2x mod 2 == 0, or, 2x mod 2 == 1."""
        DATA_OK = True
        couplings = [x1,x2,x3,x4,x5,x6]
        # Check entries for Clebsch-Gordan, Wigner 3-j, Wigner 6-j, or Racah:
        for x in couplings:
            try:
                assert ((2*x) % 2 == 0) or ((2*x) % 2 == 1)
            except AssertionError:
                logger.exception(f"Data must be entered as integral or half-integral values only: {x} is not an acceptable argument")
                DATA_OK = False
        
        if len(args) == 3:
            # Check additional data entries for Wigner 9-j:
            for y in args:
                try:
                    assert ((2*y) % 2 == 0) or ((2*y) % 2 == 1)
                except AssertionError:
                    logger.exception(f"Data must be entered as integral or half-integral values only: {y} is not an acceptable argument")
                    DATA_OK = False

        return DATA_OK

    def m_projection_rules(self):
        """Ensure magnetic-substate projections satisfy the following 
        conditions:

        (i) If (2*ji)%2 == 0, then (2*mi)%2 == 0 
            (integral ji, integral mi)

            If (2*ji)%2 == 1, then (2*mi)%2 == 1
            (half-integral ji, half-integral mi)

        (ii) |mi| <= ji

        (iii) m1 + m2 = m

        Returns:
            A `True` boolean object gets returned only if all 3 above 
            conditions are met."""
        j1, j2, j = self.j1, self.j2, self.j
        m1, m2, m = self.m1, self.m2, self.m
        
        RULE_PROJ_INT = False
        RULE_PROJ_LIMIT = False
        RULE_PROJ_SUM = False

        # Check on rule (i):
        num_passes_r1 = 0
        if (2*j1) % 2 == 0:
            try:
                assert (2*m1) % 2 == 0
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j1={self.j1} is integral; m1={self.m1} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1) requires half-integral `mi` (2*mi mod 2 = 1).")
        elif (2*j1) % 2 == 1:
            try:
                assert (2*m1) % 2 == 1
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j1={self.j1} is half-integral; m1={self.m1} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1)requires half-integral `mi` (2*mi mod 2 = 1).")

        if (2*j2) % 2 == 0:
            try:
                assert (2*m2) % 2 == 0
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j2={self.j2} is integral; m2={self.m2} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1)requires half-integral `mi` (2*mi mod 2 = 1).")
        elif (2*j2) % 2 == 1:
            try:
                assert (2*m2) % 2 == 1
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j2={self.j2} is half-integral; m2={self.m2} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1)requires half-integral `mi` (2*mi mod 2 = 1).")

        if (2*j) % 2 == 0:
            try:
                assert (2*m) % 2 == 0
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j={self.j} is integral; m={self.m} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1)requires half-integral `mi` (2*mi mod 2 = 1).")
        elif (2*j) % 2 == 1:
            try:
                assert (2*m) % 2 == 1
                num_passes_r1 += 1
            except AssertionError:
                logger.exception(f"j={self.j} is half-integral; m={self.m} is not.")
                logger.warning("If `ji` is integral (2*ji mod 2 = 0) then so must be its projection `mi` (2*mi mod 2 = 0), likewise half-integral `ji` (2*ji mod 2 = 1)requires half-integral `mi` (2*mi mod 2 = 1).")
                
        if num_passes_r1 == 3:
            RULE_PROJ_INT = True

        # Check on rule (ii):
        num_passes_r2 = 0
        try:
            assert m1 <= j1
            num_passes_r2 += 1
        except AssertionError:
            logger.exception(f"m1={self.m1} is not less than or equal to j1={self.j1}")

        try:
            assert m2 <= j2
            num_passes_r2 += 1
        except AssertionError:
            logger.exception(f"m2={self.m2} is not less than or equal to j2={self.j2}")

        try:
            assert m <= j
            num_passes_r2 += 1
        except AssertionError:
            logger.exception(f"m={self.m} is not less than or equal to j={self.j}")

        if num_passes_r2 == 3:
            RULE_PROJ_LIMIT = True

        # Check on rule (iii):
        try:
            # Test to find if "m1 + m2 = m"
            assert np.fabs((m1 + m2) - m) < np.finfo(np.float32).eps
            RULE_PROJ_SUM = True
        except AssertionError: 
            logger.exception("m1 + m2 != m")

        if RULE_PROJ_INT == True and RULE_PROJ_LIMIT == True and RULE_PROJ_SUM == True:
            return True
        else:
            logger.warning("Conditions for coupling magnetic substate projections not satisfied.")
            return

    def triangle_rule(self):
        """The triangle inequality theorem is analogous to the geometric 
        triangle inequality rule and defines the allowed range given by the sum 
        and absolute difference of angular momentum vectors j1 and j2 in 
        coupling to form total angular momentum j:
        
        |j1 - j2| <= j <= j1 + j2

        This inequality ensures that the sum of any two angular momenta is 
        greater than or equal to the third.

        Returns:
            A `True` boolean object gets returned only if the triangle 
            inequality theorem is satisfied."""
        j1, j2, j = self.j1, self.j2, self.j

        j1j2_diff = False
        j1j2_sum = False

        try:
            assert abs(j1 - j2) <= j
            j1j2_diff = True
        except AssertionError:
            logger.exception(f"Triangle inequality rule violated: |j1={self.j1} - j2={self.j2}| is not less than or equal to total angular momentum j={self.j}")
            logger.warning("Angular momentum coupling schemes must satisfy triangle inequality theorem: |j1 - j2| <= j <= j1 + j2")

        try:
            assert j <= j1 + j2
            j1j2_sum = True
        except AssertionError:
            logger.exception(f"Triangle inequality rule violated: j1={self.j1} + j2={self.j2} is not greater than or equal to total angular momentum j={self.j}")
            logger.warning("Angular momentum coupling schemes must satisfy triangle inequality theorem: |j1 - j2| <= j <= j1 + j2")

        if j1j2_diff == True and j1j2_sum == True:
            return True
        else:
            logger.warning("Conditions for coupling angular momentum vectors not satisfied.")
            return
        
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
        
        numerator = Factorial.factorial_n((j1+j2)-j) * Factorial.factorial_n((j1-j2)+j) * Factorial.factorial_n((-j1)+j2+j)
        denominator = Factorial.factorial_n(j1+j2+j+1)
        delta = sqrt(numerator/denominator)
        #print("delta_j = {0}".format(delta))
        return delta
        
    def run_v(self):
        """Function to return a tuple corresponding to the summation limits 
        needed for the angular momentum coupling scheme."""
        j1, j2, j = self.j1, self.j2, self.j
        m1, m2 = self.m1, self.m2
    
        start_0 = 0
        #start_1 = -(j-j2+m1)
        #start_2 = -(j-j1-m2)
        start_1 = j2-j-m1
        start_2 = j1-j+m2
    
        stop_0 = j1+j2-j
        stop_1 = j1-m1
        stop_2 = j2+m2
    
        #start = max(start_0, max(start_1,start_2))
        #stop = min(stop_0, min(stop_1,stop_2))
        start = max(start_0, start_1,start_2)
        stop = min(stop_0, stop_1,stop_2)
    
        return (start, stop)

    
    def cg_calc(self):
        """Evaluate Clebsch-Gordan coefficient <j1 m1 j2 m2 | j m>
        Steps involved:

        (i)   Check data input arguments for j and m.
        (ii)  Check m satisfy selection rules.
        (iii) Ensure triangle rule is satisfied for vector coupling.
        (iv)  Calculate phase factor.
        (v)   Calculate coefficient using floating point arithmetic
              (small j) and accumulators or decimal arithmetic fallback
              for large j.
        (vi)  Note that in both the primary or fallback loop the summation
              limits are determined by calling a separate function.
        """
        j1, m1 = self.j1, self.m1
        j2, m2 = self.j2, self.m2
        j, m = self.j, self.m

        data_ok = ClebschGordan.data_check(j1, m1, j2, m2, j, m)
        if data_ok is not True:
            logger.error("Input arguments must be given as integral or half-integral values only.")
            return

        m_pass = ClebschGordan.m_projection_rules(self)
        j_pass = ClebschGordan.triangle_rule(self)
        if m_pass is not True or j_pass is not True:
            logger.info("Angular momentum coupling rules not satisfied.")
            return

        dm = ClebschGordan.delta_m(self)
        if dm == 0:
            return 0.0

        # First attempt: pure float path for small j
        try:
            nom = Factorial.factorial_n(j1+m1) \
                * Factorial.factorial_n(j1-m1) \
                * Factorial.factorial_n(j2+m2) \
                * Factorial.factorial_n(j2-m2) \
                * Factorial.factorial_n(j+m) \
                * Factorial.factorial_n(j-m) \
                * ((2*j)+1) \
                * Factorial.factorial_n((j1+j2)-j) \
                * Factorial.factorial_n((j1-j2)+j) \
                * Factorial.factorial_n((-j1)+j2+j)
            den = Factorial.factorial_n(j1+j2+j+1)
            pre_factor = sqrt(nom/den)

            start = ClebschGordan.run_v(self)[0]
            stop = ClebschGordan.run_v(self)[1]
            sum_accumulator = 0.0
            for v in range(int(start), int(stop)+1):
                try:
                    dA = Factorial.factorial_n(v) \
                        * Factorial.factorial_n(((j1+j2)-j)-v) \
                        * Factorial.factorial_n((j1-m1)-v) \
                        * Factorial.factorial_n((j2+m2)-v)
                    dB = Factorial.factorial_n((j-j2)+m1+v) \
                        * Factorial.factorial_n(((j-j1)-m2)+v)
                    sum_accumulator += ((-1)**v) / (dA * dB)
                except ValueError:
                    logging.debug("break loop as soon as an argument in a factorial becomes a negative value")
                    logging.info(f"v = {v} generates negative argument in factorial")
                    break

            return float(dm * pre_factor * sum_accumulator)

        except OverflowError:
            logging.debug("Fallback into Decimal path for large-value vector coupling.")
            pass

        # Fallback: unified Decimal path
        getcontext().prec = 1000
        nom = Decimal(Factorial.factorial_n(j1+m1)) \
            * Decimal(Factorial.factorial_n(j1-m1)) \
            * Decimal(Factorial.factorial_n(j2+m2)) \
            * Decimal(Factorial.factorial_n(j2-m2)) \
            * Decimal(Factorial.factorial_n(j+m)) \
            * Decimal(Factorial.factorial_n(j-m)) \
            * Decimal((2*j)+1) \
            * Decimal(Factorial.factorial_n((j1+j2)-j)) \
            * Decimal(Factorial.factorial_n((j1-j2)+j)) \
            * Decimal(Factorial.factorial_n((-j1)+j2+j))
        den = Decimal(Factorial.factorial_n(j1+j2+j+1))
        pre_factor = (nom/den).sqrt()

        start = ClebschGordan.run_v(self)[0]
        stop = ClebschGordan.run_v(self)[1]
        sum_accumulator = Decimal(0)
        for v in range(int(start), int(stop)+1):
            try:
                dA = Decimal(Factorial.factorial_n(v)) \
                    * Decimal(Factorial.factorial_n(((j1+j2)-j)-v)) \
                    * Decimal(Factorial.factorial_n((j1-m1)-v)) \
                    * Decimal(Factorial.factorial_n((j2+m2)-v))
                dB = Decimal(Factorial.factorial_n((j-j2)+m1+v)) \
                    * Decimal(Factorial.factorial_n(((j-j1)-m2)+v))
                sum_accumulator += Decimal((-1)**v) / (dA * dB)
            except ValueError:
                logging.debug("break loop as soon as an argument in a factorial becomes a negative value")
                logging.info(f"v = {v} generates negative argument in factorial")
                break

        return float(Decimal(dm) * pre_factor * sum_accumulator)          
    
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

        data_ok = ClebschGordan.data_check(self.j1,self.j2,self.j,
                                           self.m1,self.m2,self.m)
        if data_ok == True:
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
        else:
            logger.error("Input arguments must be given as integral or half-integral values only.")

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
        
    def triangle_rule(j1,j2,j3):
        """The triangle inequality theorem is analogous to the geometric 
        triangle inequality rule and defines the allowed range given by the sum 
        and absolute difference of angular momentum vectors j1 and j2 in 
        coupling to form total angular momentum j3:
        
        |j1 - j2| <= j3 <= j1 + j2

        This inequality ensures that the sum of any two angular momenta is 
        greater than or equal to the third.

        Returns:
            A `True` boolean object gets returned only if the triangle 
            inequality theorem is satisfied."""

        j1j2_diff = False
        j1j2_sum = False

        try:
            assert abs(j1 - j2) <= j3
            j1j2_diff = True
        except AssertionError:
            logger.exception(f"Triangle inequality rule violated: |j1={j1} - j2={j2}| is not less than or equal to total angular momentum j3={j3}")
            logger.warning("Angular momentum coupling schemes must satisfy triangle inequality theorem: |j1 - j2| <= j3 <= j1 + j2")

        try:
            assert j3 <= j1 + j2
            j1j2_sum = True
        except AssertionError:
            logger.exception(f"Triangle inequality rule violated: j1={j1} + j2={j2} is not greater than or equal to total angular momentum j3={j3}")
            logger.warning("Angular momentum coupling schemes must satisfy triangle inequality theorem: |j1 - j2| <= j3 <= j1 + j2")

        if j1j2_diff == True and j1j2_sum == True:
            return True
        else:
            logger.warning("Conditions for coupling angular momentum vectors not satisfied.")
            return
        
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
        
        numerator = Factorial.factorial_n((a+b)-c) \
            * Factorial.factorial_n((a-b)+c) \
            * Factorial.factorial_n((-a)+b+c)
        denominator = Factorial.factorial_n(int(a+b+c+1))
        delta = sqrt(numerator/denominator)
        return delta

    def tri_factor_exp(a, b, c):
        """Triangular delta using log-factorials to avoid overflow.
        Returns exponentiated logarithm of the triangular delta factor."""
        # log of the numerator factorials
        log_num = (lgamma((a+b)-c+1) + lgamma((a-b)+c+1) 
                   + lgamma((-a)+b+c+1))
        # log of the denominator factorial
        log_den = lgamma(a+b+c+2)   # lgamma(n+1) == log(n!)
        # result is sqrt(exp(log_num - log_den))
        return exp(0.5 * (log_num - log_den))

    def tri_factor_log(a, b, c):
        """Return log of the triangular delta factor, to defer exponentiation."""
        log_num = (lgamma((a+b)-c+1) + lgamma((a-b)+c+1)
                   + lgamma((-a)+b+c+1))
        log_den = lgamma(a+b+c+2)
        return 0.5 * (log_num - log_den)
    
    def racah_calc(self):
        """Routine to calculate Racah or Wigner-6j without breaking up the
        calculation.  Performs `tri_factor_product` * `w` in a unified routine.

        `tri_factor_product`: This part refers to the evaluation of the product
                              of four triangular factors.
        `w`: This part part encodes the accumulator in the algebraic summation."""

        # First calculate `delta_product` routine
        triad_pass_j1j2j3 = Racah.triangle_rule(self.j1,self.j2,self.j3)
        triad_pass_j5j4j3 = Racah.triangle_rule(self.j5,self.j4,self.j3)
        triad_pass_j1j5j6 = Racah.triangle_rule(self.j1,self.j5,self.j6)
        triad_pass_j2j4j6 = Racah.triangle_rule(self.j2,self.j4,self.j6)

        if triad_pass_j1j2j3 == True and triad_pass_j5j4j3 == True and triad_pass_j1j5j6 == True and triad_pass_j2j4j6 == True:

            #delta_j1j2j3 = Racah.tri_factor(self.j1, self.j2, self.j3)
            delta_j1j2j3 = Decimal(Factorial.factorial_n((self.j1+self.j2)-self.j3)) \
                * Decimal(Factorial.factorial_n((self.j1-self.j2)+self.j3)) \
                * Decimal(Factorial.factorial_n((-self.j1)+self.j2+self.j3)) \
                / Decimal(Factorial.factorial_n(int(self.j1+self.j2+self.j3+1)))
            #delta = sqrt(numerator/denominator)
            
            #delta_j5j4j3 = Racah.tri_factor(self.j5, self.j4, self.j3)
            delta_j5j4j3 = Decimal(Factorial.factorial_n((self.j5+self.j4)-self.j3)) \
                * Decimal(Factorial.factorial_n((self.j5-self.j4)+self.j3)) \
                * Decimal(Factorial.factorial_n((-self.j5)+self.j4+self.j3)) \
                / Decimal(Factorial.factorial_n(int(self.j5+self.j4+self.j3+1)))
            #delta = sqrt(numerator/denominator)

            #delta_j1j5j6 = Racah.tri_factor(self.j1, self.j5, self.j6)
            delta_j1j5j6 = Decimal(Factorial.factorial_n((self.j1+self.j5)-self.j6)) \
                * Decimal(Factorial.factorial_n((self.j1-self.j5)+self.j6)) \
                * Decimal(Factorial.factorial_n((-self.j1)+self.j5+self.j6)) \
                / Decimal(Factorial.factorial_n(int(self.j1+self.j5+self.j6+1)))
            #delta = sqrt(numerator/denominator)
            
            #delta_j2j4j6 = Racah.tri_factor(self.j2, self.j4, self.j6)
            delta_j2j4j6 = Decimal(Factorial.factorial_n((self.j2+self.j4)-self.j6)) \
                * Decimal(Factorial.factorial_n((self.j2-self.j4)+self.j6)) \
                * Decimal(Factorial.factorial_n((-self.j2)+self.j4+self.j6)) \
                / Decimal(Factorial.factorial_n(int(self.j2+self.j4+self.j6+1)))
            #delta = sqrt(numerator/denominator)
            
            #delta_J = Newton.isqrt(delta_j1j2j3 * delta_j5j4j3 * delta_j1j5j6 * delta_j2j4j6)
            #delta_J = sqrt(delta_j1j2j3 * delta_j5j4j3 * delta_j1j5j6 * delta_j2j4j6)
            delta_J = (delta_j1j2j3 * delta_j5j4j3 * delta_j1j5j6 * delta_j2j4j6).sqrt()

            # Second, calculate `w` routine
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
        
            #w_coeff = 0
            w_coeff = Decimal(0)
            #OVERFLOW = False
            for z in range(z_min, z_max+1, 1):
                #try:
                    #numerator = (-1)**(z+b1) * Factorial.factorial_n(z+1)
                    #denominator = Factorial.factorial_n(z-a1)*Factorial.factorial_n(z-a2)*Factorial.factorial_n(z-a3)*Factorial.factorial_n(z-a4)*Factorial.factorial_n(b1-z)*Factorial.factorial_n(b2-z)*Factorial.factorial_n(b3-z)

                    #ratio = numerator/denominator
                    #w_coeff += ratio
                #except OverflowError:
                #OVERFLOW = True
                getcontext().prec = 1000
                numerator = Decimal((-1)**(z+b1)) * Decimal(Factorial.factorial_n(z+1))
                denominator = Decimal(Factorial.factorial_n(z-a1))*Decimal(Factorial.factorial_n(z-a2))*Decimal(Factorial.factorial_n(z-a3))*Decimal(Factorial.factorial_n(z-a4))*Decimal(Factorial.factorial_n(b1-z))*Decimal(Factorial.factorial_n(b2-z))*Decimal(Factorial.factorial_n(b3-z))
                
                ratio = numerator/denominator
                #ratio = float(ratio)
                w_coeff += ratio

            # Return `tri_factor_product` * `w`
            return float(delta_J * w_coeff)
    
    def phase(self):
        """Evaluate phase factor."""
        return (-1)**(self.j1+self.j2+self.j4+self.j5)
    
    def W(self):
        """Evaluate Racah coefficient."""
        data_ok = ClebschGordan.data_check(self.j1,self.j2,self.j3,
                                           self.j4,self.j5,self.j6)
        if data_ok is True:
            try:
                return Racah.racah_calc(self)
            except TypeError:
                logger.exception(f"Unspported multiplication of types: {type(Racah.delta_product(self))} and {type(Racah.w(self))}")
        else:
            logger.error("Input arguments must be given as integral or half-integral values only.")
                
    
    def symbol_6j(self):
        """Evaluate Wigner 6-j symbol."""
        data_ok = ClebschGordan.data_check(self.j1,self.j2,self.j3,
                                           self.j4,self.j5,self.j6)
        if data_ok is True:
            try:
                return Racah.W(self)*Racah.phase(self)
            except TypeError:
                logger.exception(f"Unspported multiplication of types: {type(Racah.W(self))} and {type(Racah.phase(self))}")
        else:
            logger.error("Input arguments must be given as integral or half-integral values only.")
            

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

    def __init__(self,j1,j2,j3,j4,j5,j6,j7,j8,j9,level):
        self.j1, self.j2, self.j3 = j1, j2, j3 
        self.j4, self.j5, self.j6 = j4, j5, j6 
        self.j7, self.j8, self.j9 = j7, j8, j9
        self.level=level

    def symbol_9j(self):
        """Evaluate Wigner 9-j symbol."""

        # Loop bounds are intentionally wider than needed; out-of-range h values
        # produce a TypeError which is caught to cleanly handle the triangular
        # selection rules at both ends of the summation.
        
        data_ok = ClebschGordan.data_check(self.j1,self.j2,self.j3,
                                           self.j4,self.j5,self.j6,
                                           self.j7,self.j8,self.j9)

        if data_ok is True:
        
            h_max = int(min(self.j1+self.j9,self.j4+self.j8,self.j2+self.j6)*2)
            h_min = h_max % 2
            sum_triple6j_prod = 0
            #print(f"min={h_min}, max={h_max}")
            #with LogLevelContext(logging.WARNING):
            success = False
            with LogLevelContext(self.level):
                for g in range(h_min, h_max+1):
                    h = h_min + 2*(g-1)
                    try:
                        W1 = Racah(self.j1,self.j2,self.j3,self.j6,self.j9,h/2)
                        W2 = Racah(self.j4,self.j5,self.j6,self.j2,h/2,self.j8)
                        W3 = Racah(self.j7,self.j8,self.j9,h/2,self.j1,self.j4)

                        sum_triple6j_prod = sum_triple6j_prod + (-1)**(h) * (h+1) * W1.symbol_6j() * W2.symbol_6j() * W3.symbol_6j()
                        
                        success = True

                    except TypeError:
                        if not sum_triple6j_prod:
                            continue
                        else:
                            break

            if success == True:
                return sum_triple6j_prod
            else:
                if self.level is logging.CRITICAL:
                    logger.error(f"Can't evaluate desired 9-j symbol.  Run again with debugging information to find specific problems:\n PyGammaRAD.symb9j({self.j1},{self.j2},{self.j3},{self.j4},{self.j5},{self.j6},{self.j7},{self.j8},{self.j9},logging.DEBUG)")
                    return
        else:
            logger.error("Input arguments must be given as integral or half-integral values only.")
