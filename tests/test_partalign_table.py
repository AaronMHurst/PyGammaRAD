import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class PartialAlignmentPartialTableTests(unittest.TestCase):

    __doc__="""This module contains unit tests for generation of DataFrames 
    containing tabulated partial-alignment anisotropy coefficients according to 
    the prescription outlined in Der Mateosian [1] for a defined list of 
    k-orders and Gaussian-width parameters at intervals of SIGMA/J=0.1 up to 
    a defined limit.

    References:
        [1]: E. Der Mateosian and A.W. Sunyar, At. Data and Nucl. Data
             Tables 13, 391-406 (1974).
    """

    # Unit-tests for `get_partial_table` method
    def test_get_partial_table_returns_DataFrame_default_properties(self):
        J = [10, 5/2]
        for j in J:
            df = am.get_partial_table(j)
            self.assertIsInstance(df, pd.core.frame.DataFrame)
            assert len(df) == 20

            ldf = df.values.tolist()
            for row in ldf:
                self.assertIsInstance(row, list)
                assert len(row) == 4
            
                spin = row[0]
                sJ = row[1]
                alpha2 = row[2]
                alpha4 = row[3]

                if int(2*spin) % 2 == 0:
                    self.assertIsInstance(spin, int)
                elif int(2*spin) % 2 == 1:
                    self.assertIsInstance(spin, float)
                self.assertIsInstance(sJ, float)
                self.assertIsInstance(alpha2, str)
                self.assertIsInstance(alpha4, str)

            ldf_tofloat = df.values.astype(float).tolist()
            for rowf in ldf_tofloat:
                self.assertIsInstance(rowf, list)
                assert len(rowf) == 4
                for valuef in rowf:
                    self.assertIsInstance(valuef, float)

    def test_get_partial_table_returns_DataFrame_user_properties(self):
        J = [111/2, 52]
        sJ_limit = 5.0
        k_orders = [0,2,4,6,8,10,12,14,16,18,20]
        for j in J:
            df = am.get_partial_table(j,sJ_limit,k_orders)
            self.assertIsInstance(df, pd.core.frame.DataFrame)
            assert len(df) == int(sJ_limit*10)

            ldf = df.values.tolist()
            for row in ldf:
                self.assertIsInstance(row, list)
                assert len(row) == len(k_orders)+2
            
                spin = row[0]
                sJ = row[1]

                if int(2*spin) % 2 == 0:
                    self.assertIsInstance(spin, int)
                elif int(2*spin) % 2 == 1:
                    self.assertIsInstance(spin, float)
                self.assertIsInstance(sJ, float)
                for alpha in row[2:]:
                    self.assertIsInstance(alpha, str)

            ldf_tofloat = df.values.astype(float).tolist()
            for rowf in ldf_tofloat:
                self.assertIsInstance(rowf, list)
                assert len(rowf) == len(k_orders)+2
                for valuef in rowf:
                    self.assertIsInstance(valuef, float)
            

    # Type tests
    
    def test_get_partial_table_raises_ValueError_with_wrong_value(self):
        # Non-integral or non-half-integral spin
        with self.assertRaises(ValueError):
            am.get_partial_table(2.3)
        with self.assertRaises(ValueError):
            am.get_partial_table(1.7,1,[2])

    def test_get_partial_table_raises_TypeError_with_wrong_type(self):
        # Non-numerical spin
        with self.assertRaises(TypeError):
            am.get_partial_table('x')
        with self.assertRaises(TypeError):
            am.get_partial_table('x',1.9,[2,4])
        # Non-numerical Gaussian width
        with self.assertRaises(TypeError):
            am.get_partial_table(2,'x',[2,4,6])
        with self.assertRaises(TypeError):
            am.get_partial_table(1.5,'x',[2,4,6,8])
        # Non-list object for k-orders
        with self.assertRaises(TypeError):
            am.get_partial_table(2,3.2,2)
        with self.assertRaises(TypeError):
            am.get_partial_table(2,3.2,2.0)
        with self.assertRaises(TypeError):
            am.get_partial_table(2,3.2,'2')
            
