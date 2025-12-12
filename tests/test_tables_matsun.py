import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class DerMateosianSunyarTableTests(unittest.TestCase):

    __doc__="""This module contains unit tests for generating of the tables of 
    partial-alignment anisotropy coefficients from Tables I and II in the 
    Der Mateosian paper [1].

    References:
        [1]: E. Der Mateosian and A.W. Sunyar, At. Data and Nucl. Data
             Tables 13, 391-406 (1974).
    """

    # Unit-tests for `get_table_alpha<i>` methods
    def test_get_table_alpha1_integral_J(self):
        # Retrieve data and verify attributes from integral-J table
        df = am.get_table_alpha1()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 520
        
        ldf = df.values.tolist()
        for row in ldf:
            self.assertIsInstance(row, list)
            assert len(row) == 4

            J = row[0]
            self.assertIsInstance(J, int)
            assert int(2*J)%2 == 0

            sJ = row[1]
            self.assertIsInstance(sJ, float)

            alpha2 = row[2]
            self.assertIsInstance(alpha2, str)

            alpha4 = row[3]
            self.assertIsInstance(alpha4, str)

        ldf_float = df.values.astype(float).tolist()
        for rowf in ldf_float:
            for value in rowf:
                self.assertIsInstance(value, float)

    def test_get_table_alpha2_halfint_J(self):
        # Retrieve data and verify attributes from half-integral-J table
        df = am.get_table_alpha2()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 500
        
        ldf = df.values.tolist()
        for row in ldf:
            self.assertIsInstance(row, list)
            assert len(row) == 4

            J = row[0]
            self.assertIsInstance(J, float)
            assert int(2*J)%2 == 1

            sJ = row[1]
            self.assertIsInstance(sJ, float)

            alpha2 = row[2]
            self.assertIsInstance(alpha2, str)

            alpha4 = row[3]
            self.assertIsInstance(alpha4, str)

        ldf_float = df.values.astype(float).tolist()
        for rowf in ldf_float:
            for value in rowf:
                self.assertIsInstance(value, float)
    
    # Type tests

    def test_get_table_alpha1_raises_TypeError_passing_args(self):
        with self.assertRaises(TypeError):
            am.get_table_alpha1(100)
        with self.assertRaises(TypeError):
            am.get_table_alpha1("100")
        with self.assertRaises(TypeError):
            am.get_table_alpha1([100])
        with self.assertRaises(TypeError):
            am.get_table_alpha1({"NASDAQ":100})

    def test_get_table_alpha2_raises_TypeError_passing_args(self):
        with self.assertRaises(TypeError):
            am.get_table_alpha2(500)
        with self.assertRaises(TypeError):
            am.get_table_alpha2("500")
        with self.assertRaises(TypeError):
            am.get_table_alpha2([500])
        with self.assertRaises(TypeError):
            am.get_table_alpha2({"S&P":500})

    
