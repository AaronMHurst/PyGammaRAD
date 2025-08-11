import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class ClebschGordanTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the ClebschGordan class for 
    the evaluation of Clebsch-Gordan coefficients resulting from the coupling 
    of angular momentum vectors.

    Several tests have been written to confirm well-known evaluated 
    Clebsh-Gordan results published in Table 1 of Ref. [Stevenson 2002].

    [Stevenson 2002] P.D. Stevenson, Comp. Phys. Comm. 147, 853 (2002).
    """

    # Tests confirming known evaluated Clebsch-Gordan coefficients
    
    def test_cg_confirm_value1(self):
        cgc = am.cg(2,1,1,-1,1,0)
        expected_cgc = (1/2)*np.sqrt(6/5)
        assert cgc == pytest.approx(expected_cgc)

    def test_cg_confirm_value2(self):
        cgc = am.cg(2.5,1.5,2.5,-0.5,1,1)
        expected_cgc = -2*np.sqrt(2/35)
        assert cgc == pytest.approx(expected_cgc)

    def test_cg_confirm_value3(self):
        cgc = am.cg(1,0,5,0,5,0)
        expected_cgc = 0
        assert cgc == pytest.approx(expected_cgc)

    def test_cg_confirm_value4(self):
        cgc = am.cg(3,2,2.5,-1.5,1.5,0.5)
        expected_cgc = -1*np.sqrt(1/21)
        assert cgc == pytest.approx(expected_cgc)

    def test_cg_confirm_value5(self):
        cgc = am.cg(2.5,0.5,2,0,2.5,0.5)
        expected_cgc = -2*np.sqrt(2/35)
        assert cgc == pytest.approx(expected_cgc)

    def test_cg_confirm_value6(self):
        cgc = am.cg(2.5,0.5,2,-1,2.5,-0.5)
        expected_cgc = 0
        assert cgc == pytest.approx(expected_cgc)

    # Type tests

    def test_cg_returns_float_type(self):
        cgc = am.cg(2.5,1.5,2.5,-0.5,1,1)
        self.assertIsInstance(cgc, float)

    def test_cg_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.cg()
        # Too few args
        with self.assertRaises(TypeError):
            am.cg(2,1,1)
        # Too many args
        with self.assertRaises(TypeError):
            am.cg(2,1,1,-1,1,0,-2)

    def test_cg_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.cg("2","1","1","-1","1","0")

    
