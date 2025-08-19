import pytest
import unittest
import numpy as np
import pandas as pd
from math import sqrt, factorial, gamma

import PyGammaRAD as pg
am = pg.AngularMomentum()

class FactorialTests(unittest.TestCase):

    __doc__="""This module contains unit tests to validate and compare the 
    results of factorial calculations using different methods: native, 
    recursive, and generalized.
    """
    def test_math_factorial_is_getting_called_correctly_from_class(self):
        # Test from 0 to 10,000
        for f in range(0,10000):
            assert factorial(f) == pg.Factorial.factorial_n(f)

    def test_gamma_factorial_is_getting_called_correctly_from_class(self):
        # Test from 1 to 170
        for f in range(1,170):
            assert gamma(f+1) == pg.Factorial.factorial_gamma(f)

    def test_math_factorial_returns_integer(self):
        # Test from 0 to 10,000
        for f in range(0,10000):
            fact = pg.Factorial.factorial_n(f)
            self.assertIsInstance(fact, int)

    def test_recursive_factorial_returns_integer(self):
        # Test from 0 to 2950
        for f in range(0,1000):
            fact = pg.Factorial.factorial_recursive(f)
            self.assertIsInstance(fact, int)

    def test_gamma_factorial_returns_float(self):
        # Test from 1 to 170
        for f in range(1,170):
            fact = pg.Factorial.factorial_gamma(float(f))
            self.assertIsInstance(fact, float)

    def test_math_factorial_equivalent_to_recursive_factorial(self):
        # Test from 0 to 1000
        for f in range(0,1000):
            fact_math = pg.Factorial.factorial_n(f)
            fact_recursive = pg.Factorial.factorial_recursive(f)
            assert fact_math == fact_recursive

    def test_math_factorial_equivalent_to_gamma_factorial(self):
        # Test from 1 to 23
        for f in range(1,23):
            fact_math = pg.Factorial.factorial_n(f)
            fact_gamma = pg.Factorial.factorial_gamma(float(f))
            assert fact_math == int(fact_gamma)

    def test_math_factorial_approx_same_as_gamma_factorial(self):
        # Test from 23 to 170
        ratio_expected = 1.000000000000000
        for f in range(23,170):
            fact_math = pg.Factorial.factorial_n(f)
            fact_gamma = pg.Factorial.factorial_gamma(float(f))
            assert (fact_math/fact_gamma) == pytest.approx(ratio_expected, abs=1e-15)

    def test_recursive_factorial_equivalent_to_gamma_factorial(self):
        # Test from 1 to 23
        for f in range(1,23):
            fact_recursive = pg.Factorial.factorial_recursive(float(f))
            fact_gamma = pg.Factorial.factorial_gamma(float(f))
            assert fact_gamma == fact_recursive

    def test_recursive_factorial_approx_same_as_gamma_factorial(self):
        # Test from 23 to 170
        ratio_expected = 1.000000000000000
        for f in range(23,170):
            fact_recursive = pg.Factorial.factorial_n(float(f))
            fact_gamma = pg.Factorial.factorial_gamma(float(f))
            assert (fact_recursive/fact_gamma) == pytest.approx(ratio_expected, abs=1e-15)
            
    
    
