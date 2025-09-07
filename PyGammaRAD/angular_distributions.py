from .tables import *
from .am_formulae import *
from .am_methods import *

class AngularDistributions(AngularMomentumCalculations):    
    __doc__="""Class containing functions from Yamazaki's paper [1] to 
    calculate coefficients for gamma-ray angular distributions assuming 
    complete nuclear alignment.

    References:
        [1]: T. Yamazaki, Nucl. Data Sect. A, Vol. 3, Num. 1 (1967)."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def calc_B(self, k, J):
        """Calculate the statistical population tensor B assuming complete 
        nuclear alignment for given value of k and J in accordance with Eq. (6) 
        from Yamazaki's paper [1].

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            J: A number object (int or float) representing the spin.

        Returns:
            A float value corresponding to Bk(J) which can be compared to the 
            corresponding value from Table 1 in Yamazaki's paper [1].

        Raises:
            Fewer than 2 positional arguments ('k', 'J') raises a TypeError 
            exception.

        Example:
            To calculate B(k=2,J=4):
            > calc_B(2,4)
        """
        self.k = k
        self.J = J

        if int(self.k) % 2 == 0:
            if (self.J % 2 == 0) or (self.J % 2 == 1):
                # Integral (m=0) spin
                CG = ClebschGordan(self.J, 0, self.J, 0, self.k, 0)
                cgc = CG.cg_calc()

                spin_factor = ((-1)**self.J) * np.sqrt((2*self.J) + 1)

                try:
                    B = spin_factor * cgc
                    return B
                except TypeError:
                    logger.exception(f"Unsupported multiplication of types: {type(spin_factor)} and {type(cgc)}\nReturn value is 0.0")
                    return 0.0

            elif ((float(self.J) % 2 > 1.0) and (float(self.J) % 2 < 2.0)) or ((float(self.J) % 2 > 0.0) and (float(self.J) % 2 < 1.0)):
                # Half-integral (m=0.5) spin
                CG = ClebschGordan(self.J, 0.5, self.J, -0.5, self.k, 0)
                cgc = CG.cg_calc()

                spin_factor = ((-1)**(self.J-0.5)) * np.sqrt((2*self.J) + 1)

                try:
                    B = spin_factor * cgc
                    return B
                except TypeError:
                    logger.exception(f"Unsupported multiplication of types: {type(spin_factor)} and {type(cgc)}\nReturn value is 0.0")
                    return 0.0
        else:
            logger.warning("k must be integral and even: k>0")
            return

    def calc_F(self, k, Jf, L1, L2, Ji):
        """Calculate F distribution coefficient for given value of k, Jf, L1, 
        L2, and Ji in accordance with Eq. (4) from Yamazaki's paper [1].

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            Jf: A number object (int or float) representing the final spin.
            L1: An integer object representing the order of the first 
                multipole.
            L2: An integer object representing the order of the second 
                multipole.
            Ji: A number object (int or float) representing the initial spin.
        
        Returns:
            A float value corresponding to Fk(Jf L1 L2 Ji) which can be 
            compared to the corresponding value from Table 2 in Yamazaki's 
            paper [1].

        Raises:
            Fewer than 5 positional arguments ('k', 'Jf', 'L1', 'L2', 'Ji') 
            raises a TypeError exception.

        Example:
            To calculate F(k=2,Jf=2,L1=2,L2=2,Ji=0):
            > calc_F(2,0,2,2,2)
        """
        self.k = k
        self.Jf = Jf
        self.L1 = L1
        self.L2 = L2
        self.Ji = Ji

        if int(self.k) % 2 == 0:
            
            parity_factor = (-1)**(self.Jf - self.Ji - 1)
            spin_factor = np.sqrt(((2*self.L1)+1) * ((2*self.L2)+1) * ((2*self.Ji)+1))

            CG = ClebschGordan(self.L1, 1, self.L2, -1, self.k, 0)
            cgc = CG.cg_calc()
            RC = Racah(self.Ji, self.Ji, self.k, self.L2, self.L1, self.Jf)
            W = RC.W()
            #print("Fk: W(",self.k,self.Jf,self.L1,self.L2,self.Ji,")=",W)

            try:
                F = parity_factor * spin_factor * cgc * W
                return F
            except TypeError:
                logger.exception(f"Unsupported multiplication of types: {type(parity_factor)} and {type(spin_factor)} and {type(cgc)} and {type(W)}\nReturn value is 0.0")
                return 0.0
        
        else:
            logger.warning("k must be integral and even: k>0")
            return

    def calc_BF(self, k, Jf, L1, L2, Ji):
        """Calculate BF coefficient for given value of k, Jf, L1, L2, and Ji in 
        accordance with Eq. (8), f=BF, from Yamazaki's paper [1].

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            Jf: A number object (int or float) representing the final spin.
            L1: An integer object representing the order of the first 
                multipole.
            L2: An integer object representing the order of the second 
                multipole.
            Ji: A number object (int or float) representing the initial spin.
        
        Returns:
            A float value corresponding to BF = Bk(Ji)Fk(Jf L1 L2 Ji) which 
            can be compared to the corresponding value from Table 2 in 
            Yamazaki's paper [1].

        Raises:
            Fewer than 5 positional arguments ('k', 'Jf', 'L1', 'L2', 'Ji') 
            raises a TypeError exception.

        Example:
            To calculate BF(k=2,Jf=2,L1=2,L2=2,Ji=0):
            > calc_BF(2,0,2,2,2)
        """
        self.k = k
        self.Jf = Jf
        self.L1 = L1
        self.L2 = L2
        self.Ji = Ji

        AM = AngularDistributions()
        B = AM.calc_B(self.k, self.Ji)
        F = AM.calc_F(self.k, self.Jf, self.L1, self.L2, self.Ji)
        BF = B * F
        return BF

    def calc_u(self, k, Ji, L1, Jf):
        """Calculate u coefficient uk(Ji L1 Jf) for given value of k, Ji, L1, 
        and Ji in accordance with Eq. (14) from Yamazaki's paper [1].

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            Ji: A number object (int or float) representing the initial spin.
            L1: An integer object representing the order of the unique 
                multipole.
            Jf: A number object (int or float) representing the final spin.
        
        Returns:
            A float value corresponding to uk(Ji L1 Jf) which can be compared 
            to the corresponding value from Table 2 in Yamazaki's paper [1].

        Raises:
            Fewer than 4 positional arguments ('k', 'Ji', 'L1', 'Jf') 
            raises a TypeError exception.

        Example:
            To calculate u(k=2,Ji=2,L1=2,Jf=4):
            > calc_u(2,2,2,4)
        """
        self.k = k
        self.Ji = Ji
        self.L1 = L1
        self.Jf = Jf

        if int(self.k) % 2 == 0:        
            parity_factor = (-1)**(self.Ji + self.Jf - self.L1)
            spin_factor = np.sqrt(((2*self.Ji)+1) * ((2*self.Jf)+1))
            RC = Racah(self.Ji, self.Ji, self.k, self.Jf, self.Jf, self.L1)
            W = RC.W()
            #print("uk: W(",self.k,self.Ji,self.L1,self.Jf,")=",W)

            try:
                u = parity_factor * spin_factor * W
                return u
            except TypeError:
                logger.exception(f"Unsupported multiplication of types: {type(parity_factor)} and {type(spin_factor)} and {type(W)}\nReturn value is 0.0")
                return 0.0
        else:
            logger.warning("k must be integral and even: k>0")
            return

    def A_max(self, k, Ji, L1, L2, Jf, dg=0.0):
        """Calculate the anisotropy coefficient A assuming maximum nuclear 
        alignment cf. Eq. (7) from Yamazaki's paper [1].  This calculation 
        assumes the statistical population tensor is that given by Eq. (6) for 
        complete alignment.

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            Ji: A number object (int or float) representing the initial spin.
            L1: An integer object representing the order of the first 
                multipole.
            L2: An integer object representing the order of the second 
                multipole.
            Jf: A number object (int or float) representing the final spin.
            dg: Gamma-ray mixing ratio; by default a value dg = 0.0 is assumed.
        
        Returns:
            A float value corresponding to max. Ak(Ji L1 L2 Jf).

        Raises:
            Fewer than 5 positional arguments ('k', 'Ji', 'L1', 'L2', 'Jf') 
            raises a TypeError exception.

        Example:
            To calculate A(k=2,Ji=2,L1=2,L2=2,Jf=0,dg=0):
            > A_max(2,2,2,2,0)
        """
        self.k = k
        self.Ji = Ji
        self.L1 = L1
        self.L2 = L2
        self.Jf = Jf
        self.dg = dg

        Ak = None
        AM = AngularDistributions()
        if int(self.k) % 2 == 0:
            BkFk_L1L1 = AM.calc_BF(self.k, self.Jf, self.L1, self.L1, self.Ji)
            BkFk_L1L2 = AM.calc_BF(self.k, self.Jf, self.L1, self.L2, self.Ji)
            BkFk_L2L2 = AM.calc_BF(self.k, self.Jf, self.L2, self.L2, self.Ji)
            Ak = (1/(1+(dg**2))) * (BkFk_L1L1 + (2*dg*BkFk_L1L2) + (dg**2*BkFk_L2L2))

            if (self.L1==self.L2) and (self.dg>0):
                logger.warning("Careful: Pure multipole transition defined with mixing ratio dg > 0")
                logger.warning("Result may not be physical.")
            
            if not self.dg and (self.L1!=self.L2):
                logger.info("Careful: Mixed-multipole transition defined with mixing ratio dg = 0")
                logger.info("Result reduces to B{0}({1})F{2}({3} {4} {5} {6}))".
                            format(self.k, self.Ji, self.k, self.Jf,
                                   self.L1, self.L1, self.Ji))
        
            return Ak
        else:
            logger.warning("k must be integral and even: k>0")
            return

    def U_coeff(self, k, Ji, L1, L2, Jf, dg=0.0):
        """Calculate the angular distribution coefficient Uk(Ji L1 L2 Jf) 
        cf. Eq. (13) from Yamazaki's paper [1].

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            k: An integer object representing the order.
            Ji: A number object (int or float) representing the initial spin.
            L1: An integer object representing the order of the first 
                multipole.
            L2: An integer object representing the order of the second 
                multipole.
            Jf: A number object (int or float) representing the final spin.
            dg: Gamma-ray mixing ratio; by default a value dg = 0.0 is assumed.
        
        Returns:
            A float value corresponding to Uk(Ji L1 L2 Jf).

        Raises:
            Fewer than 5 positional arguments ('k', 'Ji', 'L1', 'L2', 'Jf') 
            raises a TypeError exception.

        Example:
            To calculate U(k=2,Ji=4,L1=2,L2=2,Jf=4,dg=0):
            > U_coeff(2,4,2,2,2)
        """
        self.k = k
        self.Ji = Ji
        self.L1 = L1
        self.L2 = L2
        self.Jf = Jf
        self.dg = dg

        Ak = None
        AM = AngularDistributions()

        if int(self.k) % 2 == 0:
            uk_JiL1Jf = AM.calc_u(self.k, self.Ji, self.L1, self.Jf)
            uk_JiL2Jf = AM.calc_u(self.k, self.Ji, self.L2, self.Jf)

            Uk = (1/(1+(dg**2))) * (uk_JiL1Jf + ((dg**2)*uk_JiL2Jf))
            
            if not self.dg:
                logger.info("No mixing ratio given dg=0")
                logger.info("Result reduces to u{0}({1} {2} {3})".
                            format(self.k, self.Ji, self.L1, self.Jf))
        
            return Uk

        else:
            logger.warning("k must be integral and even: k>0")
            return
    
    def dist_W(self,A2,A4=0,*args):
        """Calculate angular distribution according to function given by 
        Eq. (2) in Yamazaki's paper [1].  Either theoretical (Ak) or 
        experimental (ak) anisotropy coefficients can be passed as arguments.

        Notes:
            [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).

        Arguments:
            A2: An float object correspodning to the k=2 anisotropy 
                coefficient.
            A4: An float object correspodning to the k=4 anisotropy 
                coefficient.
            args: An optional numpy array object corresponding to a user-defined
                  angular range.

        Returns:
            A tuple object of two array-type elements:

            [0]: Calculated angular distribution (numpy.ndarray);
            [1]: Angular range over which angular distribution is calculated
                 (numpy.ndarray).

        Raises:
            Passing a non-numpy.ndarray object as an <*args> parameter 
            raises an AssertionError exception.

            Failure to pass 2 positional arguments for 'A2' and 'A4' 
            raises a TypeError exception.

        Examples:
            To calculate W(theta) over default range (0,180,1800) assuming 
            maximum alignment for 2+->0+ stretched E2 quadrupole transition:
            > dist_W(A_max(2,2,2,2,0),A_max(4,2,2,2,0))

            To calculate the same W(theta) above assuming maximum alignment 
            over a user-defined range <numpy.ndarray>:
            > dist_W(A_max(2,2,2,2,0),A_max(4,2,2,2,0),<numpy.ndarray>)

            To calculate W(theta) over defualt range for a2=0.47, a4=-0.41:
            > dist_W(0.47,-0.41)
        """

        self.A2 = A2
        self.A4 = A4
        theta = None
        P2, P4 = None, None

        L = Legendre()
        if len(args) == 0 or args == []:
            theta = L.deg2rad()
            P2 = L.lpoly2()
            P4 = L.lpoly4()
        elif len(args) == 1:
            theta = args[0]
            try:
                assert type(theta) is np.ndarray
                P2 = L.lpoly2(theta)
                P4 = L.lpoly4(theta)
            except AssertionError:
                logger.exception("The angular range must be passed as a 1D numpy array\n{0} is not a numpy.ndarray object".format(type(theta)))
                return
        else:
            logger.error("Too many arguments", exc_info=False)
            logger.warning("The angular range should be passed as a 1D numpy array")
            return

        W = 1 + (A2*P2) + (A4*P4)

        return W, theta
    
class Legendre(AngularDistributions):
    __doc__="""Legendre polynomials: Pk as a function of cos(theta)."""

    def __init__(self):
        self.theta = np.linspace(0.0,180.0,1800)

    def deg2rad(self,x=None):
        """Method enabling the conversion of arrays, lists and single-valued 
        floats and integers from degrees to radians.

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.

        Returns:
            The corresponding angle or set of angles in radians as a numpy 
            array or floating-point value depending on the input arguments.

        Examples:
            (i) To convert the default set of angles from 0 to 180 degrees to 
            radians within the initailized linear space (0,180,1800):
            > am.deg2rad()

            (ii) To convert a numpy array <numpy.ndarray> or <list> of values 
            from degrees to radians
            > am.deg2rad(<numpy.ndarray>)
            > am.deg2rad(<list>)
            
            (iii) To convert a single-value float <float> or <integer> from 
            degress to radians:
            > am.deg2rad(<float>)
            > am.deg2rad(<int>)
        """
        if x is None:
            return (self.theta/180.0)*np.pi
        else:
            if type(x) is list:
                x = np.array(x)
            return (x/180.0)*np.pi

    def rad2deg(self,x=None):
        """Method enabling the conversion of arrays, lists and single-valued 
        floats and integers from radians to degrees.

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.

        Returns:
            The corresponding angle or set of angles in degrees as a numpy 
            array or floating-point value depending on the input arguments.

        Examples:
            (i) To return the default set of angles from 0 to 180 degrees 
            within the initailized linear space (0,180,1800):
            > am.deg2rad()

            (ii) To convert a numpy array <numpy.ndarray> or <list> of values 
            from radians to degrees
            > am.rad2deg(<numpy.ndarray>)
            > am.rad2deg(<list>)
            
            (iii) To convert a single-value float <float> or <integer> from 
            radians to degrees:
            > am.rad2deg(<float>)
            > am.rad2deg(<int>)
        """
        if x is None:
            #L = Legendre()
            #return L.rad2deg(L.deg2rad())
            return self.theta
        else:
            if type(x) is list:
                x = np.array(x)
            return (180.0/np.pi)*x
        
    def lpoly0(self,x=None):
        """Legendre polynomial of degree k=0 (zeroth order) as a function of 
        the cosine of the angle (theta): P0(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=0 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=0 Legendre polynomal:
            > lpoly0()
 
            (i) To get the k=0 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly0(<numpy.ndarray>)
            > lpoly0(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly0(90)
        """
        # Suppress RuntimeWarninig to avoid NaN results from zero division in
        # array elements
        np.seterr(invalid='ignore')
        if x is None:
            #return self.theta/self.theta
            # Only divide non-zeroes; replace NaN from zero division with 1.0
            return np.where(self.theta!=0,np.divide(self.theta,self.theta),1.0)
        else:
            #return x/x
            if (type(x) is float) or (type(x) is int):
                if not x:
                    return 1.0
                else:
                    return x/x
            else:
                return np.where(x!=0, np.divide(x,x), 1.0)

    def lpoly1(self,x=None):
        """Legendre polynomial of degree k=1 (first order) as a function of 
        the cosine of the angle (theta): P1(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=1 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=1 Legendre polynomal:
            > lpoly1()
 
            (i) To get the k=1 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly1(<numpy.ndarray>)
            > lpoly1(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly1(90)
        """
        L = Legendre()
        if x is None:
            return np.cos(L.deg2rad())
        else:
            if type(x) is list:
                x = np.array(x)
            return np.cos(L.deg2rad(x))
        
    def lpoly2(self,x=None):
        """Legendre polynomial of degree k=2 (second order) as a function of 
        the cosine of the angle (theta): P2(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=2 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=2 Legendre polynomal:
            > lpoly2()
 
            (i) To get the k=2 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly2(<numpy.ndarray>)
            > lpoly2(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly2(90)
        """
        L = Legendre()
        if x is None:
            return 0.5 * (3*np.cos(L.deg2rad())**2 - 1)
        else:
            if type(x) is list:
                x = np.array(x)
            return 0.5 * (3*np.cos(L.deg2rad(x))**2 - 1)

    def lpoly3(self,x=None):
        """Legendre polynomial of degree k=3 (third order) as a function of 
        the cosine of the angle (theta): P3(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=3 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=3 Legendre polynomal:
            > lpoly3()
 
            (i) To get the k=3 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly3(<numpy.ndarray>)
            > lpoly3(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly3(90)
        """
        L = Legendre()
        if x is None:
            return 0.5 * ((5*np.cos(L.deg2rad())**3) - (3*np.cos(L.deg2rad())))
        else:
            if type(x) is list:
                x = np.array(x)
            return 0.5 * ((5*np.cos(L.deg2rad(x))**3) - (3*np.cos(L.deg2rad(x))))

    def lpoly4(self,x=None):
        """Legendre polynomial of degree k=4 (fourth order) as a function of 
        the cosine of the angle (theta): P4(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=4 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=4 Legendre polynomal:
            > lpoly4()
 
            (i) To get the k=4 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly4(<numpy.ndarray>)
            > lpoly4(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly4(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/8.0) * ((35*np.cos(L.deg2rad())**4) - (30*np.cos(L.deg2rad())**2) + 3)
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/8.0) * ((35*np.cos(L.deg2rad(x))**4) - (30*np.cos(L.deg2rad(x))**2) + 3)

    def lpoly5(self,x=None):
        """Legendre polynomial of degree k=5 (fifth order) as a function of 
        the cosine of the angle (theta): P5(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=5 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=5 Legendre polynomal:
            > lpoly5()
 
            (i) To get the k=5 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly5(<numpy.ndarray>)
            > lpoly5(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly5(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/8.0) * ((63*np.cos(L.deg2rad())**5) - (70*np.cos(L.deg2rad())**3) + (15*np.cos(L.deg2rad())))
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/8.0) * ((63*np.cos(L.deg2rad(x))**5) - (70*np.cos(L.deg2rad(x))**3) + (15*np.cos(L.deg2rad(x))))

    def lpoly6(self,x=None):
        """Legendre polynomial of degree k=6 (sixth order) as a function of 
        the cosine of the angle (theta): P6(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=6 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=6 Legendre polynomal:
            > lpoly6()
 
            (i) To get the k=6 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly6(<numpy.ndarray>)
            > lpoly6(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly6(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/16.0) * ((231*np.cos(L.deg2rad())**6) - (315*np.cos(L.deg2rad())**4) + (105*np.cos(L.deg2rad())**2) - 5)
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/16.0) * ((231*np.cos(L.deg2rad(x))**6) - (315*np.cos(L.deg2rad(x))**4) + (105*np.cos(L.deg2rad(x))**2) - 5)

    def lpoly7(self,x=None):
        """Legendre polynomial of degree k=7 (seventh order) as a function of 
        the cosine of the angle (theta): P7(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=7 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=7 Legendre polynomal:
            > lpoly7()
 
            (i) To get the k=7 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly7(<numpy.ndarray>)
            > lpoly7(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly7(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/16.0) * ((429*np.cos(L.deg2rad())**7) - (693*np.cos(L.deg2rad())**5) + (315*np.cos(L.deg2rad())**3) - (35*np.cos(L.deg2rad())))
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/16.0) * ((429*np.cos(L.deg2rad(x))**7) - (693*np.cos(L.deg2rad(x))**5) + (315*np.cos(L.deg2rad(x))**3) - (35*np.cos(L.deg2rad(x))))
            
    def lpoly8(self,x=None):
        """Legendre polynomial of degree k=8 (eighth order) as a function of 
        the cosine of the angle (theta): P8(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=8 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=8 Legendre polynomal:
            > lpoly8()
 
            (i) To get the k=8 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly8(<numpy.ndarray>)
            > lpoly8(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly8(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/128.0) * ((6435*np.cos(L.deg2rad())**8) - (12012*np.cos(L.deg2rad())**6) + (6930*np.cos(L.deg2rad())**4) - (1260*np.cos(L.deg2rad())**2) + 35)
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/128.0) * ((6435*np.cos(L.deg2rad(x))**8) - (12012*np.cos(L.deg2rad(x))**6) + (6930*np.cos(L.deg2rad(x))**4) - (1260*np.cos(L.deg2rad(x))**2) + 35)

    def lpoly9(self,x=None):
        """Legendre polynomial of degree k=9 (ninth order) as a function of 
        the cosine of the angle (theta): P9(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=9 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=9 Legendre polynomal:
            > lpoly9()
 
            (i) To get the k=9 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly9(<numpy.ndarray>)
            > lpoly9(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly9(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/128.0) * ((12155*np.cos(L.deg2rad())**9) - (25740*np.cos(L.deg2rad())**7) + (18018*np.cos(L.deg2rad())**5) - (4620*np.cos(L.deg2rad())**3) + (315*np.cos(L.deg2rad())))
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/128.0) * ((12155*np.cos(L.deg2rad(x))**9) - (25740*np.cos(L.deg2rad(x))**7) + (18018*np.cos(L.deg2rad(x))**5) - (4620*np.cos(L.deg2rad(x))**3) + (315*np.cos(L.deg2rad(x))))

    def lpoly10(self,x=None):
        """Legendre polynomial of degree k=10 (tenth order) as a function of 
        the cosine of the angle (theta): P10(cos(theta)).

        Arguments:
            None: Default linear space initialized will be used: 0,180,1800 
            x: User-defined linear space may also be passed as a 1D numpy array
               or list object. Single-valued numerical float or integer objects 
               are also acceptable arguments.  All values entered in units of 
               degrees.

        Returns:
            The legendre polynomial of degree k=3 over the desired angular 
            (theta)range as either an numpy array object containing 
            floating-point object elements, or single-valued float depending on
            input arguments.
        
        Examples:
            (i) To get the default k=10 Legendre polynomal:
            > lpoly10()
 
            (i) To get the k=10 Legendre polynomal for a user-defined range 
            <numpy.ndarray> or <list>:
            > lpoly10(<numpy.ndarray>)
            > lpoly10(<list>)

            (iii) To get the value of the Legendre polynomial at 90 degrees:
            > lpoly10(90)
        """
        L = Legendre()
        if x is None:
            return (1.0/256.0) * ((46189*np.cos(L.deg2rad())**10) - (109395*np.cos(L.deg2rad())**8) + (90090*np.cos(L.deg2rad())**6) - (30030*np.cos(L.deg2rad())**4) + (3465*np.cos(L.deg2rad())**2) - 63)
        else:
            if type(x) is list:
                x = np.array(x)
            return (1.0/256.0) * ((46189*np.cos(L.deg2rad(x))**10) - (109395*np.cos(L.deg2rad(x))**8) + (90090*np.cos(L.deg2rad(x))**6) - (30030*np.cos(L.deg2rad(x))**4) + (3465*np.cos(L.deg2rad(x))**2) - 63)

    
