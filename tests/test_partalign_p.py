import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class PartialAlignmentStatisticalPkTests(unittest.TestCase):

    __doc__="""This module contains unit tests for the statistical tensor 
    function according to Eq. (1) by Yamazaki [1] or Eq (4) by 
    Der Mateosian [2].

    References:
        [1]: T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
        [2]: E. Der Mateosian and A.W. Sunyar, At. Data and Nucl. Data
             Tables 13, 391-406 (1974).
    """
    def test_partial_p_results_and_ensure_float_returns_integral_J(self):
        # Test reproduction of tensors for range of k-orders for integral J
        # Test to ensure return of floating-point value type
        p_k_2_to_20_J_10 = [-0.8481100972980309, 0.46578818760623686,
                            -0.19065618723392358, 0.05979929687660883,
                            -0.014494167334627597, 0.0026921531520451613,
                            -0.00037233023238878365, 3.604830876877075e-05,
                            -2.1333754664957916e-06, 5.0240332875500975e-08]
        
        k_orders = [k for k in range(2,21,2)]
        J = 10
        sJ = 0.3
        for i, k in enumerate(k_orders):
            p = am.partial_p(k, J, sJ)
            self.assertIsInstance(p, float)
            for j, p_known in enumerate(p_k_2_to_20_J_10):
                if i == j:
                    assert p == pytest.approx(p_known, abs=1.0e-10)
    
    def test_partial_p_results_and_ensure_float_returns_halfint_J(self):
        # Test reproduction of tensors for range of k-orders for half-integral J
        # Test to ensure return of floating-point value type
        p_k_2_to_20_J_75_2 = [-0.8268540124180678, 0.43717675931830724,
                              -0.1743470492887113, 0.055261974470383786,
                              -0.014430888972202844, 0.00318671461958084,
                              -0.0006067854042506346, 0.00010111679226054882,
                              -1.4918226963003459e-05, 1.966316329057277e-06]

        k_orders = [k for k in range(2,21,2)]
        J = 75/2
        sJ = 0.3
        for i, k in enumerate(k_orders):
            p = am.partial_p(k, J, sJ)
            self.assertIsInstance(p, float)
            for j, p_known in enumerate(p_k_2_to_20_J_75_2):
                if i == j:
                    assert p == pytest.approx(p_known, abs=1.0e-10)
        
        
    # Type tests
    
    def test_partial_p_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.partial_p()
        # Too few args
        with self.assertRaises(TypeError):
            am.partial_p(2)
        with self.assertRaises(TypeError):
            am.partial_p(2,10)
        # Too many args
        with self.assertRaises(TypeError):
            am.partial_p(2,10,0.3,2)
        with self.assertRaises(TypeError):
            am.partial_p(2,10,0.3,2,7.5)

    def test_partial_p_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.partial_p("2","10","0.3")

           
    def test_partial_p_raises_ValueError_with_k_odd_arg(self):
        # Pass odd-k args
        with self.assertRaises(ValueError):
            am.partial_p(1,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_p(3,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_p(5,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_p(7,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_p(9,10,0.3)
        with self.assertRaises(ValueError):
            am.partial_p(11,10,0.3)

