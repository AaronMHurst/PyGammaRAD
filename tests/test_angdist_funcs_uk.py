import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionUkTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the angular 
    distribution function for the Uk-coefficient calculated in accordance with 
    Equation (14) from Yamazaki's paper [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_calc_u_equal_to_table2a_values_and_returns_float(self):
        # Test to compare calculated U coefficients with tabulated values.
        # Loop over integral spins from Table 2(a)
        # Test to ensure return of floating-point value type
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_orders = [2, 4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                U2_tabulated = row[8]
                U4_tabulated = row[9]

                try:
                    if int(L1) == int(L2):
                        Uk_calculated = am.calc_u(k, Ji, L1, Jf)
                        if k == 2:
                            assert Uk_calculated == pytest.approx(U2_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_calculated, float)
                
                        elif k == 4:
                            assert Uk_calculated == pytest.approx(U4_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_calculated, float)
                except ValueError:
                    if int(U2_tabulated) == 0 or int(U4_tabulated) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_u(k, Ji, L1, Jf)

    
    def test_calc_u_equal_to_table2b_values_and_returns_float(self):
        # Test to compare calculated U coefficients with tabulated values.
        # Loop over half-integral spins from Table 2(b)
        # Test to ensure return of floating-point value type
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_orders = [2, 4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                U2_tabulated = row[8]
                U4_tabulated = row[9]

                try:
                    if int(L1) == int(L2):
                        Uk_calculated = am.calc_u(k, Ji, L1, Jf)
                        if k == 2:
                            assert Uk_calculated == pytest.approx(U2_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_calculated, float)
                        elif k == 4:
                            assert Uk_calculated == pytest.approx(U4_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_calculated, float)
                except ValueError:
                    if int(U2_tabulated) == 0 or int(U4_tabulated) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_u(k, Ji, L1, Jf)


    # Test Yamazaki vs Rose & Brink methods for uk(Ji L1 Jf) determination

    def test_calc_u_Yamazaki_equal_to_Racah_ratio_RoseBrink_integral_J(self):
        # Test to compare calculated U coefficients using Eq. 14 [Yamazaki]
        # gives same result as Racah ratio defined in Eq. 3.45 [Rose and Brink].
        # Loop over integral spins from Table 2(a) to get input arguments for
        # equations.
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_orders = [2, 4, 6, 8, 10]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                try:
                    if int(L1) == int(L2):
                        # Calculate uk using Yamazaki's method
                        Uk_yamazaki = am.calc_u(k, Ji, L1, Jf)

                        # Calculate uk from Rose & Brinks Racah ratios
                        Uk_rosebrink_1 = (-1)**(k) * am.racah(Ji, Jf, Ji, Jf, L1, k)/am.racah(Ji, Jf, Ji, Jf, L1, 0)
                        Uk_rosebrink_2 = (-1)**(k) * am.racah(Ji, Ji, Jf, Jf, k, L1)/am.racah(Ji, Ji, Jf, Jf, 0, L1)

                        assert Uk_yamazaki == pytest.approx(Uk_rosebrink_1, abs=1.0e-09)
                        assert Uk_yamazaki == pytest.approx(Uk_rosebrink_2, abs=1.0e-09)
                        assert Uk_rosebrink_1 == pytest.approx(Uk_rosebrink_2, abs=1.0e-09)
                        if int(Ji)==2 and int(Jf)==4 and int(L1)==2:
                            # Specific test case
                            u2 = 0.74915
                            u4 = 0.28472
                            if k == 2:
                                assert Uk_yamazaki == pytest.approx(u2, abs=1.0e-05)
                                assert Uk_rosebrink_1 == pytest.approx(u2, abs=1.0e-05)
                                assert Uk_rosebrink_2 == pytest.approx(u2, abs=1.0e-05)
                            if k == 4:
                                assert Uk_yamazaki == pytest.approx(u4, abs=1.0e-05)
                                assert Uk_rosebrink_1 == pytest.approx(u4, abs=1.0e-05)
                                assert Uk_rosebrink_2 == pytest.approx(u4, abs=1.0e-05)

                except TypeError:
                    with self.assertRaises(TypeError):
                        am.racah(Ji, Jf, Ji, Jf, L1, k)/am.racah(Ji, Jf, Ji, Jf, L1, 0)

    def test_calc_u_Yamazaki_equal_to_Racah_ratio_RoseBrink_halfint_J(self):
        # Test to compare calculated U coefficients using Eq. 14 [Yamazaki]
        # gives same result as Racah ratio defined in Eq. 3.45 [Rose and Brink].
        # Loop over half-integral spins from Table 2(b) to get input arguments
        # for equations.
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_orders = [2, 4, 6, 8, 10]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                try:
                    if int(L1) == int(L2):
                        # Calculate uk using Yamazaki's method
                        Uk_yamazaki = am.calc_u(k, Ji, L1, Jf)

                        # Calculate uk from Rose & Brinks Racah ratios
                        Uk_rosebrink_1 = (-1)**(k) * am.racah(Ji, Jf, Ji, Jf, L1, k)/am.racah(Ji, Jf, Ji, Jf, L1, 0)
                        Uk_rosebrink_2 = (-1)**(k) * am.racah(Ji, Ji, Jf, Jf, k, L1)/am.racah(Ji, Ji, Jf, Jf, 0, L1)

                        assert Uk_yamazaki == pytest.approx(Uk_rosebrink_1, abs=1.0e-09)
                        assert Uk_yamazaki == pytest.approx(Uk_rosebrink_2, abs=1.0e-09)
                        assert Uk_rosebrink_1 == pytest.approx(Uk_rosebrink_2, abs=1.0e-09)
                        if int(2*Ji)==9 and int(2*Jf)==9 and int(L1)==1:
                            # Specific test case
                            u2 = 0.87879
                            u4 = 0.59596
                            if k == 2:
                                assert Uk_yamazaki == pytest.approx(u2, abs=1.0e-05)
                                assert Uk_rosebrink_1 == pytest.approx(u2, abs=1.0e-05)
                                assert Uk_rosebrink_2 == pytest.approx(u2, abs=1.0e-05)
                            if k == 4:
                                assert Uk_yamazaki == pytest.approx(u4, abs=1.0e-05)
                                assert Uk_rosebrink_1 == pytest.approx(u4, abs=1.0e-05)
                                assert Uk_rosebrink_2 == pytest.approx(u4, abs=1.0e-05)

                except TypeError:
                    with self.assertRaises(TypeError):
                        am.racah(Ji, Jf, Ji, Jf, L1, k)/am.racah(Ji, Jf, Ji, Jf, L1, 0)
                    
    # Type tests
    
    def test_calc_u_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_u()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_u(4,2,2)
        # Too many args
        with self.assertRaises(TypeError):
            am.calc_u(4,2,2,4,2)

    def test_calc_u_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_u("4","2","2","4")
    
    def test_calc_u_returns_NoneType_with_k_odd_arg(self):
        U = am.calc_u(1,2,2,4)
        self.assertIsNone(U)
        U = am.calc_u(3,2,2,4)
        self.assertIsNone(U)
        U = am.calc_u(5,2,2,4)
        self.assertIsNone(U)
    
