import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class Wigner3jTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the Wigner3j class for the 
    evaluation of Wigner-3j symbols resulting from the coupling of angular 
    momentum vectors.

    Several tests have been written to confirm the results of  well-known 
    evaluated 3-j symbols published in Table 2 of Ref. [Stevenson 2002].

    [Stevenson 2002] P.D. Stevenson, Comp. Phys. Comm. 147, 853 (2002).
    """

    # Tests confirming known evaluated Wigner-3j symbols

    def test_symb3j_confirm_value1(self):
        s3j = am.symb3j(1,1,1,1,0,-1)
        expected_s3j = -1/np.sqrt(6)
        assert s3j == pytest.approx(expected_s3j)

    def test_symb3j_confirm_value2(self):
        s3j = am.symb3j(3,3,3,-1,-1,2)
        expected_s3j = 0
        assert s3j == pytest.approx(expected_s3j)

    def test_symb3j_confirm_value3(self):
        s3j = am.symb3j(3.5,2.5,2,3.5,-1.5,-2)
        expected_s3j = -np.sqrt(2)/6
        assert s3j == pytest.approx(expected_s3j)

    def test_symb3j_confirm_value4(self):
        s3j = am.symb3j(7.5,7.5,0,1.5,-1.5,0)
        expected_s3j = 1/4
        assert s3j == pytest.approx(expected_s3j)

    def test_symb3j_confirm_value5(self):
        s3j = am.symb3j(8,5.5,4.5,2,-3.5,1.5)
        expected_s3j = (-127/12)*np.sqrt(6/46189)
        assert s3j == pytest.approx(expected_s3j)

    def test_symb3j_confirm_value6(self):
        s3j = am.symb3j(8,7,6,3,0,-3)
        expected_s3j = 0
        assert s3j == pytest.approx(expected_s3j)

    # Type tests
    
    def test_symb3j_returns_float_type(self):
        s3j = am.symb3j(8,5.5,4.5,2,-3.5,1.5)
        self.assertIsInstance(s3j, float)

    def test_symb3j_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.symb3j()
        # Too few args
        with self.assertRaises(TypeError):
            am.symb3j(3.5,2.5,2)
        # Too many args
        with self.assertRaises(TypeError):
            am.symb3j(3.5,2.5,2,3.5,-1.5,-2,-7.5)

    def test_symb3j_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.symb3j("1","1","1","1","0","-1")
