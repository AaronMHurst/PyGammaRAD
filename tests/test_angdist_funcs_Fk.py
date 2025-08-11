import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionFkTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the angular 
    distribution function for the Fk-coefficient calculated in accordance with 
    Equation (4) from Yamazaki's paper [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_calc_F_equal_to_table2a_values_and_returns_float(self):
        # Test to compare calculated F coefficients with tabulated values.
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

                F2_tabulated = row[4]
                F4_tabulated = row[6]

                try:
                    Fk_calculated = am.calc_F(k, Jf, L1, L2, Ji)
                    if k == 2:
                        assert Fk_calculated == pytest.approx(F2_tabulated, abs=0.00001)
                        self.assertIsInstance(Fk_calculated, float)
                    elif k == 4:
                        assert Fk_calculated == pytest.approx(F4_tabulated, abs=0.00001)
                        self.assertIsInstance(Fk_calculated, float)
                except ValueError:
                    if int(F2_tabulated) == 0 or int(F4_tabulated) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_F(k, Jf, L1, L2, Ji)

    def test_calc_F_equal_to_table2b_values_and_returns_float(self):
        # Test to compare calculated F coefficients with tabulated values.
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

                F2_tabulated = row[4]
                F4_tabulated = row[6]

                try:
                    Fk_calculated = am.calc_F(k, Jf, L1, L2, Ji)
                    if k == 2:
                        assert Fk_calculated == pytest.approx(F2_tabulated, abs=0.00001)
                        self.assertIsInstance(Fk_calculated, float)
                    elif k == 4:
                        assert Fk_calculated == pytest.approx(F4_tabulated, abs=0.00001)
                        self.assertIsInstance(Fk_calculated, float)
                except ValueError:
                    if int(F2_tabulated) == 0 or int(F4_tabulated) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_F(k, Jf, L1, L2, Ji)
    
    # Type tests
    
    def test_calc_F_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_F()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_F(2,0,2,2)
        # Too many args
        with self.assertRaises(TypeError):
            am.calc_F(2,0,2,2,2,2)

    def test_calc_F_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_F("2","0","2","2","2")
    
    def test_calc_F_returns_NoneType_with_k_odd_arg(self):
        F = am.calc_F(1,0,2,2,2)
        self.assertIsNone(F)
        F = am.calc_F(3,0,2,2,2)
        self.assertIsNone(F)
        F = am.calc_F(5,0,2,2,2)
        self.assertIsNone(F)
    
