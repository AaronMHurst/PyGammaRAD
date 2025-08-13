import pytest
import unittest
import numpy as np
import pandas as pd
from math import sqrt, factorial

import PyGammaRAD as pg
am = pg.AngularMomentum()

class NewtonIsqrtTests(unittest.TestCase):

    __doc__="""This module contains unit tests to validate the `isqrt` function 
    based on Newton's iterative method in evaluating the square root of large 
    integers.
    """

    def test_isqrt_returns_integers(self):
        for i in range(0,1000000):
            isqrt = pg.Newton.isqrt(i)
            self.assertIsInstance(isqrt, int)
        
    def test_isqrt_returns_same_result_as_math_sqrt(self):
        # Verify sqrts from 0 -> 1,000,000: math.sqrt cf. isqrt
        for i in range(0,1000000):
            math_sqrt = sqrt(i)
            self.assertIsInstance(math_sqrt, float)
            isqrt = pg.Newton.isqrt(i)
            assert isqrt == int(math_sqrt)

    def test_isqrt_returns_same_result_as_numpy_sqrt(self):
        # Verify sqrts from 0 -> 1,000,000: np.sqrt cf. isqrt
        for i in range(0,1000000):
            np_sqrt = np.sqrt(i)
            self.assertIsInstance(np_sqrt, float)
            isqrt = pg.Newton.isqrt(i)
            assert isqrt == int(np_sqrt)

    def test_isqrt_factorials_same_result_as_math_sqrt_factorials(self):
        # Test all factorials up to handling limit of math.sqrt(170!)
        sig_diff = [32,36,37,49,62,65,76,77,78,80,83,
                    85,124,130,144,146,152,158,163]
        
        for f in range(0,171):
            sqrt_fact = sqrt(factorial(f))
            isqrt_fact = pg.Newton.isqrt(factorial(f))
            diff = abs(sqrt_fact - isqrt_fact)
            if f < 32:
                assert diff < 1.0
            elif f in sig_diff:
                assert diff > 1.0
            else:
                #assert int(diff) == 0
                assert not diff

    def test_isqrt_factorials_same_result_as_numpy_sqrt_factorials(self):
        # Test all factorials up to handling limit of np.sqrt(20!)
        for f in range(0,20):
            npsqrt_fact = np.sqrt(factorial(f))
            isqrt_fact = pg.Newton.isqrt(factorial(f))
            diff = abs(npsqrt_fact - isqrt_fact)
            # Assert floating-point differences of < 1
            self.assertIsInstance(diff, float)
            assert diff < 1.0
            # Assert no integral difference
            #assert int(diff) == 0
            assert not int(diff)
            
    def test_beyond_math_method_limits_raises_OverflowError(self):
        for i in range(0,200):
            if i <= 170:
                math_sqrt = sqrt(factorial(i))
                self.assertIsInstance(math_sqrt, float)
            else:
                with self.assertRaises(OverflowError):
                    sqrt(factorial(i))

    def test_beyond_numpy_method_limits_raises_TypeError(self):
        for i in range(0,30):
            if i <= 20:
                np_sqrt = np.sqrt(factorial(i))
                self.assertIsInstance(np_sqrt, float)
            else:
                with self.assertRaises(TypeError):
                    np.sqrt(factorial(i))

    def test_isqrt_returns_int_for_very_large_factorials(self):
        # Returns integers for factorials well beyond limits of
        # traditional square root methods
        for i in range(0,1000):
            isqrt_fact = pg.Newton.isqrt(factorial(i))
            self.assertIsInstance(isqrt_fact, int)

