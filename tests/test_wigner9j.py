import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class Wigner9jTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the Wigner6j class for the 
    evaluation of Wigner-9j symbols resulting from the coupling of a set of 
    angular momenta.

    Several tests have been written to confirm the results of well-known 
    evaluated 9-j symbols published in Table 4 of Ref. [Stevenson 2002].

    Results were also checked and confirmed against Anthony Stone's Wigner 
    coefficient calculator online [Stone].

    [Stevenson 2002] P.D. Stevenson, Comp. Phys. Comm. 147, 853 (2002).
    [Stone] https://www-stone.ch.cam.ac.uk/wigner.shtml
    """

    # Tests confirming known evaluated Wigner-9j symbols

    def test_symb9j_confirm_value1(self):
        s9j = am.symb9j(1,3,2,1.5,2.5,2,0.5,0.5,0)
        #expected_s9j = (1/6)*np.sqrt(1/5) # CPC result appears to be wrong
        expected_s9j = -(1/6)*np.sqrt(1/35)
        assert s9j == pytest.approx(expected_s9j)

    def test_symb9j_confirm_value2(self):
        s9j = am.symb9j(3,4,2,3.5,3.5,2,0.5,0.5,1)
        expected_s9j = (1/3)*np.sqrt(1/210)
        assert s9j == pytest.approx(expected_s9j)

    def test_symb9j_confirm_value3(self):
        s9j = am.symb9j(1,3,4,0.5,3.5,3,0.5,0.5,1)
        expected_s9j = -1/168
        assert s9j == pytest.approx(expected_s9j)

    def test_symb9j_confirm_value4(self):
        s9j = am.symb9j(3,3.5,3.5,2.5,3,3.5,0.5,0.5,1)
        expected_s9j = (-11/336)*np.sqrt(1/231)
        assert s9j == pytest.approx(expected_s9j)

    def test_symb9j_confirm_value5(self):
        """Triads featuring maximal coupling with large angular momenta"""
        s9j = am.symb9j(20,20,40,20,20,40,20,20,40)
        expected_s9j = (6525173305508/1666331383246797314439) \
            *np.sqrt(4843218585190677539/5959445981) # [Stevenson 2002]
        #expected_s9j = (283703187196/1666331383246797314439) \
        #    *np.sqrt(111394027459385583397/259106347) 
        assert s9j == pytest.approx(expected_s9j)

    # Type tests

    def test_symb9j_returns_float_type(self):
        s9j = am.symb9j(3,4,2,3.5,3.5,2,0.5,0.5,1)
        self.assertIsInstance(s9j, float)

    def test_symb9j_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.symb9j()
        # Too few args
        with self.assertRaises(TypeError):
            am.symb9j(3,4,2,3.5,3.5,2)
        # Too many args
        with self.assertRaises(TypeError):
            am.symb9j(3,4,2,3.5,3.5,2,0.5,0.5,1,1.5)

    #def test_raises_TypeError_with_wrong_types_args(self):
    # Others raise TypeError...
    def test_symb9j_raises_ValueError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(ValueError):
            am.symb9j("3","4","2","3.5","3.5","2","0.5","0.5","1")
    
