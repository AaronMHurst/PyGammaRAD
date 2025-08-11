import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class Wigner6jTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the Wigner6j class for the 
    evaluation of Wigner-6j symbols resulting from the coupling of a set of 
    angular momenta.

    Several tests have been written to confirm the results of well-known 
    evaluated 6-j symbols published in Table 3 of Ref. [Stevenson 2002].

    [Stevenson 2002] P.D. Stevenson, Comp. Phys. Comm. 147, 853 (2002).
    """

    # Tests confirming known evaluated Wigner-6j symbols

    def test_symb6j_confirm_value1(self):
        s6j = am.symb6j(2,2,1,2,1,2)
        expected_s6j = 7/(10*np.sqrt(21))
        assert s6j == pytest.approx(expected_s6j)

    def test_symb6j_confirm_value2(self):
        s6j = am.symb6j(3.5,3,1.5,1,1.5,3)
        expected_s6j = 0
        assert s6j == pytest.approx(expected_s6j)

    def test_symb6j_confirm_value3(self):
        s6j = am.symb6j(4,4,1,4,3,1)
        expected_s6j = -1/36
        assert s6j == pytest.approx(expected_s6j)

    def test_symb6j_confirm_value4(self):
        s6j = am.symb6j(6,6,4,4.5,3.5,5.5)
        expected_s6j = 1363/36036
        assert s6j == pytest.approx(expected_s6j)

    def test_symb6j_confirm_value5(self):
        s6j = am.symb6j(6.5,6,1.5,3,3.5,6)
        expected_s6j = 0
        assert s6j == pytest.approx(expected_s6j)

    def test_symb6j_confirm_value6(self):
        s6j = am.symb6j(8,8,8,8,7,7)
        expected_s6j = (2557/193154)*np.sqrt(7/3)
        assert s6j == pytest.approx(expected_s6j)

    # Type tests
    
    def test_symb6j_returns_float_type(self):
        s6j = am.symb6j(8,8,8,8,7,7)
        self.assertIsInstance(s6j, float)

    def test_symb6j_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.symb6j()
        # Too few args
        with self.assertRaises(TypeError):
            am.symb6j(8,8,8)
        # Too many args
        with self.assertRaises(TypeError):
            am.symb6j(8,8,8,8,7,7,7)

    def test_symb6j_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.symb6j("8","8","8","8","7","7")
