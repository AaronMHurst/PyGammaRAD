import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class PartialAlignmentPopulationParameterTests(unittest.TestCase):

    __doc__="""This module contains unit tests for calculation of the 
    population parameter given by Eq. (11) in Yamazaki [1] or Eq. (6) in 
    Der Mateosian [2].

    References:
        [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
        [2]: E. Der Mateosian and A.W. Sunyar, At. Data and Nucl. Data
             Tables 13, 391-406 (1974).
    """

    # Unit-tests for `partial_P` method
    def test_partial_P_results_and_ensure_float_returns_integral_J(self):
        # Test reproduction of population parameters for all magnetic
        # substates of J=5
        # Test to ensure return of floating-point value type
        P_J_5 = [0.00102838008447911, 0.007598758135239185,
                 0.03600077212843083, 0.10936068950970002,
                 0.2130055377112537, 0.26601172486179436,
                 0.2130055377112537, 0.10936068950970002,
                 0.03600077212843083, 0.007598758135239185,
                 0.00102838008447911]

        J = 5
        sJ = 0.3
        for i, m in enumerate(range(-J,J+1,1)):
            P = am.partial_P(J, m, sJ)
            self.assertIsInstance(P, float)
            for j, m_known in enumerate(P_J_5):
                if i == j:
                    assert P == pytest.approx(m_known, abs=1.0e-10)

    def test_partial_P_results_and_ensure_float_returns_halfint_J(self):
        # Test reproduction of population parameters for all magnetic
        # substates of J=5/2
        # Test to ensure return of floating-point value type
        P_J_5_2 = [0.0020564748131369714, 0.07199155473056555,
                   0.42595197045629746, 0.42595197045629746,
                   0.07199155473056555, 0.0020564748131369714]

        J = 5/2
        sJ = 0.3
        for i, m in enumerate(range(-int(2*J),int(2*J)+1)):
            if m%2 == 1:
                m = m/2
                P = am.partial_P(J, m, sJ)
                self.assertIsInstance(P, float)
                for j, P_known in enumerate(P_J_5_2):
                    if i == (2*j):
                        assert P == pytest.approx(P_known, abs=1.0e-10)
        
    # Type tests: `partial_P`

    def test_partial_P_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.partial_P()
        # Too few args
        with self.assertRaises(TypeError):
            am.partial_P(5)
        with self.assertRaises(TypeError):
            am.partial_P(5,3)
        # Too many args
        with self.assertRaises(TypeError):
            am.partial_P(5,3,0.3,2)
        with self.assertRaises(TypeError):
            am.partial_P(5,3,0.3,2,7.5)

    def test_partial_P_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.partial_P("5","3","0.3")
    
    def test_partial_P_is_NoneType_with_m_GT_J_or_m_LT_J(self):
        # Pass m > J; integral
        P = am.partial_P(5,7,0.3)
        self.assertIsNone(P)
        # Pass m > J; half-integral 
        P = am.partial_P(5/2,7/2,0.3)
        self.assertIsNone(P)
        
        # Pass m < J; integral
        P = am.partial_P(8,-9,0.3)
        self.assertIsNone(P)
        # Pass m < J; half-integral 
        P = am.partial_P(19/2,-21/2,0.3)
        self.assertIsNone(P)

    def test_partial_P_integral_J_must_have_integral_m_projection_likewise_halfint_otherwise_AssertionError_gets_raised(self):
        # Half-integral J / Integral m
        with self.assertRaises(AssertionError):
            am.partial_P(5/2,2,0.3)
        with self.assertRaises(AssertionError):
            am.partial_P(11/2,-5,0.3)
        # Integral J / Half-integral m
        with self.assertRaises(AssertionError):
            am.partial_P(3,1/2,0.3)
        with self.assertRaises(AssertionError):
            am.partial_P(10,-7/2,0.3)

            
    # Unit-tests for `pop_paras` method
    def test_pop_paras_default_returned_properties(self):
        J = 5
        pops = am.pop_paras(J)
        # Returns tuple object containing 2 elements
        self.assertIsInstance(pops, tuple)
        assert len(pops) == 2

        # First element is list of length 2J+1
        assert len(pops[0]) == (2*J)+1
        for m in pops[0]:
            self.assertIsInstance(m, int)

        # Second element is dictionary of length 20
        assert len(pops[1]) == 20
        for sJ, P in pops[1].items():
            self.assertIsInstance(sJ, float)
            self.assertIsInstance(P, list)
            assert len(P) == (2*J)+1
            for p in P:
                self.assertIsInstance(p, float)

    def test_pop_paras_user_SIGMA_J_returned_properties(self):
        J = 21/2
        sigmaJ = [0.1,0.5,2.0,3.0,5.0]
        pops = am.pop_paras(J,sigmaJ)
        # Returns tuple object containing 2 elements
        self.assertIsInstance(pops, tuple)
        assert len(pops) == 2

        # First element is list of length 2J+1
        assert len(pops[0]) == (2*J)+1
        for m in pops[0]:
            self.assertIsInstance(m, float)

        # Second element is dictionary of length 20
        assert len(pops[1]) == len(sigmaJ)
        for sJ, P in pops[1].items():
            self.assertIsInstance(sJ, float)
            self.assertIsInstance(P, list)
            assert len(P) == (2*J)+1
            for p in P:
                self.assertIsInstance(p, float)

    # Type tests: `pop_paras`

    def test_pop_paras_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.pop_paras()
        # Too many args
        with self.assertRaises(TypeError):
            am.pop_paras(5,[3],2.0)

    def test_partial_P_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.partial_P("5","3")
        # sJ argument must be passed as a list to avoid TypeError exception
        with self.assertRaises(TypeError):
            am.pop_paras(5,3)
    
        

        
