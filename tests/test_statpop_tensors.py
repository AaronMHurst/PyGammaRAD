import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class PopulationTensorTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing calculated statistical 
    population tensors assuming complete nuclear alignment Bk(J) with those 
    obtained from Table 1 of Ref. [1967Ya05].  The tests are performed for 
    integral J and half-integral J values.  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """

    def test_calc_B_equal_to_table1_values(self):
        """Test to compare calculated B values with tabulated values."""
        # Loop over integral spins
        for J in range(1,21):
            for k in range(2,8,2):
                # Test for k=2,4,6 in nested loop
                B_table = am.get_row_table1(J)[int((k-2)/2)]
                try:
                    B_calculated = am.calc_B(k,J)
                    assert B_calculated == pytest.approx(B_table, abs=0.00001)
                except ValueError:
                    if int(B_table) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_B(k,J)

        # Loop over half-integral spins
        for J in range(1,21):
            J = J + 0.5
            for k in range(2,8,2):
                # Test for k=2,4,6 in nested loop
                B_table = am.get_row_table1(J)[int((k-2)/2)]
                try:
                    B_calculated = am.calc_B(k,J)
                    assert B_calculated == pytest.approx(B_table, abs=0.00001)
                except ValueError:
                    if int(B_table) == 0:
                        with self.assertRaises(ValueError):
                            am.calc_B(k,J)
            
    # Type tests

    def test_calc_B_returns_float_type(self):

        k = 2
        # Tensors with k=2 always return floating point values
        for J in range(1,51):
            # half-integral J
            B = am.calc_B(k,J)
            self.assertIsInstance(B, float)
                
            # half-integral J
            J += 0.5
            B = am.calc_B(k,J)
            self.assertIsInstance(B, float)

        
        for k in range(4,8,2):
            # Tensors with k=4,6 usually return floating point values
            # Exceptions: B4 J < 2; B6 J < 3
            for J in range(1,51):
                try:
                    B = am.calc_B(k,J)
                    self.assertIsInstance(B, float)
                except ValueError:
                    if (k == 4 and J < 2) or (k == 6 and J < 3):
                        with self.assertRaises(ValueError):
                            am.calc_B(k,J)
        
    def test_calc_B_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_B()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_B(2)
        # Too many args
        with self.assertRaises(TypeError):
            am.calc_B(2,10,6)

    def test_calc_B_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_B("2","10")

    def test_calc_B_returns_NoneType_with_k_odd_arg(self):
        B = am.calc_B(1,2)
        self.assertIsNone(B)
        B = am.calc_B(3,2)
        self.assertIsNone(B)
        B = am.calc_B(5,4)
        self.assertIsNone(B)
        
