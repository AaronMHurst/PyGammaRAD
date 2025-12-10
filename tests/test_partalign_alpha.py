import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class PartialAlignmentAlphaTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the partial-alignment 
    anisotropy angular distribution coefficient, `alpha`, calculated in 
    accordance with Equation (10) from Yamazaki [] to the tabulated data in 
    Der Mateosian and Sunyar [].  The calculation method for `alpha` is reliant 
    upon the correct determination of the population parameter, given by 
    Equation (10) in Yamazaki [] or Equation (6) in Der Mateosian and Sunyar [],
    and the statistical population tensor, given by Equation 1 in Yamazaki [] 
    or Equation (4) in Der Mateosian and Sunyar [], in addition to the tensor 
    for complete alignment, given by Equation (6) in Yamazaki [] or 
    Equation (5) in Der Mateosian and Sunyar [].  Accordingly, the unit-tests 
    for the `alpha` coefficient provide further verification of the calculation 
    methods underpinning its evaluation.

    [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_partial_a_equal_to_get_table_alpha1_values_and_returns_float_integral_J(self):
        # Test to compare calculated `alpha` coefficients with tabulated values.
        # Loop over integral spins from `alpha`-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_table_alpha1()
        ldf = df.values.astype(float).tolist()
        
        for row in ldf:
            J = int(row[0])
            sJ = row[1]

            alpha2_tabulated = row[2]
            alpha4_tabulated = row[3]

            try:
                alpha2_calculated = am.partial_a(2, J, sJ)
                assert alpha2_calculated == pytest.approx(alpha2_tabulated, abs=0.000001)
                self.assertIsInstance(alpha2_calculated, float)

            except TypeError:
                epsilon = 1e-09
                if alpha2_tabulated - epsilon < 0:
                    with self.assertRaises(TypeError):
                        am.partial_a(2, J, sJ)

            try:
                alpha4_calculated = am.partial_a(4, J, sJ)
                assert alpha4_calculated == pytest.approx(alpha4_tabulated, abs=0.000001)
                self.assertIsInstance(alpha4_calculated, float)

            except TypeError:
                epsilon = 1e-09
                if alpha2_tabulated - epsilon < 0:
                    with self.assertRaises(TypeError):
                        am.partial_a(4, J, sJ)

    def test_partial_a_equal_to_get_table_alpha2_values_and_returns_float_halfint_J(self):
        # Test to compare calculated `alpha` coefficients with tabulated values.
        # Loop over half-integral spins from `alpha`-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_table_alpha2()
        ldf = df.values.astype(float).tolist()
        
        for row in ldf:
            J = row[0]
            sJ = row[1]

            alpha2_tabulated = row[2]
            alpha4_tabulated = row[3]

            try:
                alpha2_calculated = am.partial_a(2, J, sJ)
                assert alpha2_calculated == pytest.approx(alpha2_tabulated, abs=0.000001)
                self.assertIsInstance(alpha2_calculated, float)

            except TypeError:
                epsilon = 1e-09
                if alpha2_tabulated - epsilon < 0:
                    with self.assertRaises(TypeError):
                        am.partial_a(2, J, sJ)

            try:
                alpha4_calculated = am.partial_a(4, J, sJ)
                assert alpha4_calculated == pytest.approx(alpha4_tabulated, abs=0.000001)
                self.assertIsInstance(alpha4_calculated, float)

            except TypeError:
                epsilon = 1e-09
                if alpha2_tabulated - epsilon < 0:
                    with self.assertRaises(TypeError):
                        am.partial_a(4, J, sJ)


    # Type tests
    
    def test_partial_a_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.partial_a()
        # Too few args
        with self.assertRaises(TypeError):
            am.partial_a(2)
        with self.assertRaises(TypeError):
            am.partial_a(2,10)
        # Too many args
        with self.assertRaises(TypeError):
            am.partial_a(2,10,0.3,2)
        with self.assertRaises(TypeError):
            am.partial_a(2,10,0.3,2,7.5)

    def test_partial_a_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.partial_a("2","10","0.3")

           
    def test_partial_a_raises_ValueError_with_k_odd_arg(self):
        # Pass odd-k args
        with self.assertRaises(ValueError):
            am.partial_a(1,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_a(3,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_a(5,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_a(7,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_a(9,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_a(11,10,0.3)
        


