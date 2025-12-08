import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionRkTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the angular 
    distribution function for the Rk-coefficient calculated in accordance with 
    Equations (3.36) and (3.37) from Rose and Brink's review article paper 

    [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_calc_R_equal_to_tableRa_values_and_returns_float_integral_J(self):
        # Test to compare calculated R coefficients with tabulated values.
        # Loop over integral spins from R-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tableRa()
        ldf = df.values.tolist()

        for row in ldf:
            Ji = row[0]
            Jf = row[1]
            L1 = row[2]
            L2 = row[3]

            R2L1L1_tabulated = row[4]
            R2L1L2_tabulated = row[5]
            R2L2L2_tabulated = row[6]

            try:
                R2L1L1_calculated = am.calc_R(2, L1, L1, Ji, Jf)
                assert R2L1L1_calculated == pytest.approx(R2L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R2L1L1_calculated, float)

            except ValueError:
                if int(R2L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L1, L1, Ji, Jf)

            try:
                R2L1L2_calculated = am.calc_R(2, L1, L2, Ji, Jf)
                assert R2L1L2_calculated == pytest.approx(R2L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R2L1L2_calculated, float)

            except ValueError:
                if int(R2L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L1, L2, Ji, Jf)

            try:
                R2L2L2_calculated = am.calc_R(2, L2, L2, Ji, Jf)
                assert R2L2L2_calculated == pytest.approx(R2L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R2L2L2_calculated, float)

            except ValueError:
                if int(R2L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L2, L2, Ji, Jf)

            R4L1L1_tabulated = row[7]
            R4L1L2_tabulated = row[8]
            R4L2L2_tabulated = row[9]

            try:
                R4L1L1_calculated = am.calc_R(4, L1, L1, Ji, Jf)
                assert R4L1L1_calculated == pytest.approx(R4L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R4L1L1_calculated, float)

            except ValueError:
                if int(R4L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L1, L1, Ji, Jf)

            try:
                R4L1L2_calculated = am.calc_R(4, L1, L2, Ji, Jf)
                assert R4L1L2_calculated == pytest.approx(R4L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R4L1L2_calculated, float)

            except ValueError:
                if int(R4L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L1, L2, Ji, Jf)

            try:
                R4L2L2_calculated = am.calc_R(4, L2, L2, Ji, Jf)
                assert R4L2L2_calculated == pytest.approx(R4L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R4L2L2_calculated, float)

            except ValueError:
                if int(R4L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L2, L2, Ji, Jf)

            R6L1L1_tabulated = row[10]
            R6L1L2_tabulated = row[11]
            R6L2L2_tabulated = row[12]

            try:
                R6L1L1_calculated = am.calc_R(6, L1, L1, Ji, Jf)
                assert R6L1L1_calculated == pytest.approx(R6L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R6L1L1_calculated, float)

            except ValueError:
                if int(R6L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L1, L1, Ji, Jf)

            try:
                R6L1L2_calculated = am.calc_R(6, L1, L2, Ji, Jf)
                assert R6L1L2_calculated == pytest.approx(R6L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R6L1L2_calculated, float)

            except ValueError:
                if int(R6L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L1, L2, Ji, Jf)

            try:
                R6L2L2_calculated = am.calc_R(6, L2, L2, Ji, Jf)
                assert R6L2L2_calculated == pytest.approx(R6L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R6L2L2_calculated, float)

            except ValueError:
                if int(R6L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L2, L2, Ji, Jf)

            R8L1L1_tabulated = row[13]
            R8L1L2_tabulated = row[14]
            R8L2L2_tabulated = row[15]

            try:
                R8L1L1_calculated = am.calc_R(8, L1, L1, Ji, Jf)
                assert R8L1L1_calculated == pytest.approx(R8L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R8L1L1_calculated, float)

            except ValueError:
                if int(R8L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L1, L1, Ji, Jf)

            try:
                R8L1L2_calculated = am.calc_R(8, L1, L2, Ji, Jf)
                assert R8L1L2_calculated == pytest.approx(R8L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R8L1L2_calculated, float)

            except ValueError:
                if int(R8L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L1, L2, Ji, Jf)

            try:
                R8L2L2_calculated = am.calc_R(8, L2, L2, Ji, Jf)
                assert R8L2L2_calculated == pytest.approx(R8L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R8L2L2_calculated, float)

            except ValueError:
                if int(R8L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L2, L2, Ji, Jf)

    def test_calc_R_equal_to_tableRb_values_and_returns_float_halfint_J(self):
        # Test to compare calculated R coefficients with tabulated values.
        # Loop over integral spins from R-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tableRb()
        ldf = df.values.tolist()

        for row in ldf:
            Ji = row[0]
            Jf = row[1]
            L1 = row[2]
            L2 = row[3]

            R2L1L1_tabulated = row[4]
            R2L1L2_tabulated = row[5]
            R2L2L2_tabulated = row[6]

            try:
                R2L1L1_calculated = am.calc_R(2, L1, L1, Ji, Jf)
                assert R2L1L1_calculated == pytest.approx(R2L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R2L1L1_calculated, float)

            except ValueError:
                if int(R2L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L1, L1, Ji, Jf)

            try:
                R2L1L2_calculated = am.calc_R(2, L1, L2, Ji, Jf)
                assert R2L1L2_calculated == pytest.approx(R2L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R2L1L2_calculated, float)

            except ValueError:
                if int(R2L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L1, L2, Ji, Jf)

            try:
                R2L2L2_calculated = am.calc_R(2, L2, L2, Ji, Jf)
                assert R2L2L2_calculated == pytest.approx(R2L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R2L2L2_calculated, float)

            except ValueError:
                if int(R2L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(2, L2, L2, Ji, Jf)

            R4L1L1_tabulated = row[7]
            R4L1L2_tabulated = row[8]
            R4L2L2_tabulated = row[9]

            try:
                R4L1L1_calculated = am.calc_R(4, L1, L1, Ji, Jf)
                assert R4L1L1_calculated == pytest.approx(R4L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R4L1L1_calculated, float)

            except ValueError:
                if int(R4L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L1, L1, Ji, Jf)

            try:
                R4L1L2_calculated = am.calc_R(4, L1, L2, Ji, Jf)
                assert R4L1L2_calculated == pytest.approx(R4L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R4L1L2_calculated, float)

            except ValueError:
                if int(R4L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L1, L2, Ji, Jf)

            try:
                R4L2L2_calculated = am.calc_R(4, L2, L2, Ji, Jf)
                assert R4L2L2_calculated == pytest.approx(R4L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R4L2L2_calculated, float)

            except ValueError:
                if int(R4L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(4, L2, L2, Ji, Jf)

            R6L1L1_tabulated = row[10]
            R6L1L2_tabulated = row[11]
            R6L2L2_tabulated = row[12]

            try:
                R6L1L1_calculated = am.calc_R(6, L1, L1, Ji, Jf)
                assert R6L1L1_calculated == pytest.approx(R6L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R6L1L1_calculated, float)

            except ValueError:
                if int(R6L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L1, L1, Ji, Jf)

            try:
                R6L1L2_calculated = am.calc_R(6, L1, L2, Ji, Jf)
                assert R6L1L2_calculated == pytest.approx(R6L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R6L1L2_calculated, float)

            except ValueError:
                if int(R6L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L1, L2, Ji, Jf)

            try:
                R6L2L2_calculated = am.calc_R(6, L2, L2, Ji, Jf)
                assert R6L2L2_calculated == pytest.approx(R6L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R6L2L2_calculated, float)

            except ValueError:
                if int(R6L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(6, L2, L2, Ji, Jf)

            R8L1L1_tabulated = row[13]
            R8L1L2_tabulated = row[14]
            R8L2L2_tabulated = row[15]

            try:
                R8L1L1_calculated = am.calc_R(8, L1, L1, Ji, Jf)
                assert R8L1L1_calculated == pytest.approx(R8L1L1_tabulated, abs=0.0001)
                self.assertIsInstance(R8L1L1_calculated, float)

            except ValueError:
                if int(R8L1L1_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L1, L1, Ji, Jf)

            try:
                R8L1L2_calculated = am.calc_R(8, L1, L2, Ji, Jf)
                assert R8L1L2_calculated == pytest.approx(R8L1L2_tabulated, abs=0.0001)
                self.assertIsInstance(R8L1L2_calculated, float)

            except ValueError:
                if int(R8L1L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L1, L2, Ji, Jf)

            try:
                R8L2L2_calculated = am.calc_R(8, L2, L2, Ji, Jf)
                assert R8L2L2_calculated == pytest.approx(R8L2L2_tabulated, abs=0.0001)
                self.assertIsInstance(R8L2L2_calculated, float)

            except ValueError:
                if int(R8L2L2_tabulated) == 0: 
                    with self.assertRaises(ValueError):
                        am.calc_R(8, L2, L2, Ji, Jf)


    
    # Type tests
    
    def test_calc_R_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_R()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_R(2,1,2,1.5)
        # Too many args
        with self.assertRaises(TypeError):
            am.calc_R(2,1,2,1.5,2.5,2.5)

    def test_calc_R_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_R("2","1","2","1.5","2.5")

    def test_calc_R_returns_NoneType_with_k_odd_arg(self):
        R = am.calc_R(1,1,2,1.5,2.5)
        self.assertIsNone(R)
        R = am.calc_R(3,1,2,1.5,2.5)
        self.assertIsNone(R)
        R = am.calc_R(5,1,2,1.5,2.5)
        self.assertIsNone(R)
    
