import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionSkTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the angular 
    distribution function for the Sk-coefficient calculated in accordance with 
    Equation (3.59) from Rose and Brink's review article paper 

    [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_calc_S_equal_to_tableSa_values_and_returns_float_integral_J(self):
        # Test to compare calculated S coefficients with tabulated values.
        # Loop over integral spins from S-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tableSa()
        ldf = df.values.tolist()

        for row in ldf:
            l1 = row[0]
            l2 = row[1]
            J = row[2]
            s = row[3]

            S0_tabulated = row[4]
            S2_tabulated = row[5]
            S4_tabulated = row[6]
            S6_tabulated = row[7]
            S8_tabulated = row[8]

            Sk_table_list = [S0_tabulated, S2_tabulated, S4_tabulated,
                             S6_tabulated, S8_tabulated]
            
            try:
                S0_calculated = am.calc_S(l1, l2, J, s, 0)
                assert S0_calculated == pytest.approx(S0_tabulated, abs=0.00001)
                self.assertIsInstance(S0_calculated, float)

            except ValueError:
                if int(S0_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 0)

            try:
                S2_calculated = am.calc_S(l1, l2, J, s, 2)
                assert S2_calculated == pytest.approx(S2_tabulated, abs=0.00001)
                self.assertIsInstance(S2_calculated, float)

            except ValueError:
                if int(S2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 2)

            try:
                S4_calculated = am.calc_S(l1, l2, J, s, 4)
                assert S4_calculated == pytest.approx(S4_tabulated, abs=0.00001)
                self.assertIsInstance(S4_calculated, float)

            except ValueError:
                if int(S4_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 4)

            try:
                S6_calculated = am.calc_S(l1, l2, J, s, 6)
                assert S6_calculated == pytest.approx(S6_tabulated, abs=0.00001)
                self.assertIsInstance(S6_calculated, float)

            except ValueError:
                if int(S6_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 6)

            try:
                S8_calculated = am.calc_S(l1, l2, J, s, 8)
                assert S8_calculated == pytest.approx(S8_tabulated, abs=0.00001)
                self.assertIsInstance(S8_calculated, float)

            except ValueError:
                if int(S8_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 8)

            # Test calling method for all k-orders
            Sk_calc_list = am.calc_S(l1, l2, J, s)
            for i,Sk_t in enumerate(Sk_table_list):
                for j,Sk_c in enumerate(Sk_calc_list):
                    if i == j:
                        assert Sk_c == pytest.approx(Sk_t, abs=0.00001)
                        self.assertIsInstance(Sk_c, float)
                        self.assertIsInstance(Sk_t, float)

    def test_calc_S_equal_to_tableSb_values_and_returns_float_halfint_J(self):
        # Test to compare calculated S coefficients with tabulated values.
        # Loop over half-integral spins from S-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tableSb()
        ldf = df.values.tolist()

        for row in ldf:
            l1 = row[0]
            l2 = row[1]
            J = row[2]
            s = row[3]

            S0_tabulated = row[4]
            S2_tabulated = row[5]
            S4_tabulated = row[6]
            S6_tabulated = row[7]
            S8_tabulated = row[8]

            Sk_table_list = [S0_tabulated, S2_tabulated, S4_tabulated,
                             S6_tabulated, S8_tabulated]
            
            try:
                S0_calculated = am.calc_S(l1, l2, J, s, 0)
                assert S0_calculated == pytest.approx(S0_tabulated, abs=0.00001)
                self.assertIsInstance(S0_calculated, float)

            except ValueError:
                if int(S0_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 0)

            try:
                S2_calculated = am.calc_S(l1, l2, J, s, 2)
                assert S2_calculated == pytest.approx(S2_tabulated, abs=0.00001)
                self.assertIsInstance(S2_calculated, float)

            except ValueError:
                if int(S2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 2)

            try:
                S4_calculated = am.calc_S(l1, l2, J, s, 4)
                assert S4_calculated == pytest.approx(S4_tabulated, abs=0.00001)
                self.assertIsInstance(S4_calculated, float)

            except ValueError:
                if int(S4_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 4)

            try:
                S6_calculated = am.calc_S(l1, l2, J, s, 6)
                assert S6_calculated == pytest.approx(S6_tabulated, abs=0.00001)
                self.assertIsInstance(S6_calculated, float)

            except ValueError:
                if int(S6_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 6)

            try:
                S8_calculated = am.calc_S(l1, l2, J, s, 8)
                assert S8_calculated == pytest.approx(S8_tabulated, abs=0.00001)
                self.assertIsInstance(S8_calculated, float)

            except ValueError:
                if int(S8_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_S(l1, l2, J, s, 8)

            # Test calling method for all k-orders
            Sk_calc_list = am.calc_S(l1, l2, J, s)
            for i,Sk_t in enumerate(Sk_table_list):
                for j,Sk_c in enumerate(Sk_calc_list):
                    if i == j:
                        assert Sk_c == pytest.approx(Sk_t, abs=0.00001)
                        self.assertIsInstance(Sk_c, float)
                        self.assertIsInstance(Sk_t, float)
                
    
    # Type tests
    
    def test_calc_S_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_S()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_S(4,4,5)

    def test_calc_S_returns_NoneType_with_too_many_args(self):
        # Too many args
        S = am.calc_S(4,4,5,4,2,2)
        self.assertIsNone(S)
        S = am.calc_S(4,4,5,4,2,2,8)
        self.assertIsNone(S)
        S = am.calc_S(4,4,5,4,2,2,8,10)
        self.assertIsNone(S)
    
    def test_calc_S_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_S("4","4","5","4")
        with self.assertRaises(TypeError):
            am.calc_S("4","4","5","4","2")

    def test_calc_S_returns_NoneType_with_k_odd_arg(self):
        S = am.calc_S(4,4,5,4,1)
        self.assertIsNone(S)
        S = am.calc_S(4,4,5,4,3)
        self.assertIsNone(S)
        S = am.calc_S(4,4,5,4,5)
        self.assertIsNone(S)
