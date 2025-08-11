import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionUkMaxTests(unittest.TestCase):
    __doc__="""This module contains unit tests to verify the calculation of the 
    Uk(Ji L1 L2 Jf) angular distribution function  assuming complete nuclear 
    alignment, in accordance with Eq. (13) from Yamazaki's paper [1].

    References:
        [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    def test_U_coeff_equivalent_tabulated_uk_integral_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(a); integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                U2_tabulated = row[8]
                U4_tabulated = row[9]

                if int(L1) == int(L2):
                    try:
                        Uk_max = am.U_coeff(k, Ji, L1, L2, Jf)
                        if k == 2:
                            assert Uk_max == pytest.approx(U2_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_max, float)
                        elif k == 4:
                            assert Uk_max == pytest.approx(U4_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_max, float)
                    except ValueError:
                        if not U2_tabulated or U4_tabulated:
                            with self.assertRaises(ValueError):
                                am.U_coeff(k, Ji, L1, L2, Jf)

    def test_U_coeff_equivalent_tabulated_uk_halfint_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(a); integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                U2_tabulated = row[8]
                U4_tabulated = row[9]

                if int(L1) == int(L2):
                    try:
                        Uk_max = am.U_coeff(k, Ji, L1, L2, Jf)
                        if k == 2:
                            assert Uk_max == pytest.approx(U2_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_max, float)
                        elif k == 4:
                            assert Uk_max == pytest.approx(U4_tabulated, abs=0.00001)
                            self.assertIsInstance(Uk_max, float)
                    except ValueError:
                        if not U2_tabulated or U4_tabulated:
                            with self.assertRaises(ValueError):
                                am.U_coeff(k, Ji, L1, L2, Jf)

    # Type tests

    def test_U_coeff_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.U_coeff()
        # Too few args
        with self.assertRaises(TypeError):
            am.U_coeff(2,4,2,2)
        # Too many args
        with self.assertRaises(TypeError):
            # default 6th argument dg=0
            am.U_coeff(2,2,2,2,2,-0.18,2)

    def test_U_coeff_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.U_coeff("2","4","2","2","2")
        with self.assertRaises(TypeError):
            am.U_coeff("2","2","2","2","2","-0.18")

    def test_U_coeff_raises_TypeError_with_k_odd_arg(self):
        # U_coeff method tries to return summation of uk coefficients
        # Odd-k attempts <float> * <NoneType>
        with self.assertRaises(TypeError):
            am.U_coeff(1,1.5,1,1,2.5)
        with self.assertRaises(TypeError):
            am.U_coeff(3,1.5,2,2,1.5)
        with self.assertRaises(TypeError):
            am.U_coeff(5,2.5,3,3,4.5)

    
