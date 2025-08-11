import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionAkMaxTests(unittest.TestCase):

    __doc__="""This module contains unit tests to compare calculated anisotropy 
    coefficients based on complete nuclear alignment to BkFk coefficients 
    listed in Tables 2(a) and 2(b) from Yamazaki's paper [1].  

    References:
    [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """

    def test_A_max_equivalent_tabulated_BF_J_integral_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(a); integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                B2F2_tabulated = row[5]
                B4F4_tabulated = row[7]

                if int(L1) == int(L2):

                    try:
                        Ak_max = am.A_max(k, Ji, L1, L2, Jf)
                        if k == 2:
                            assert Ak_max == pytest.approx(B2F2_tabulated, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                        elif k == 4:
                            assert Ak_max == pytest.approx(B4F4_tabulated, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                    except ValueError:
                        if not B2F2_tabulated or B4F4_tabulated:
                            with self.assertRaises(ValueError):
                                am.A_max(k, Ji, L1, L2, Jf)

    def test_A_max_equivalent_tabulated_BF_J_halfint_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(b); half-integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                B2F2_tabulated = row[5]
                B4F4_tabulated = row[7]

                if int(L1) == int(L2):

                    try:
                        Ak_max = am.A_max(k, Ji, L1, L2, Jf)
                        if k == 2:
                            assert Ak_max == pytest.approx(B2F2_tabulated, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                        elif k == 4:
                            assert Ak_max == pytest.approx(B4F4_tabulated, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                    except ValueError:
                        if not B2F2_tabulated or B4F4_tabulated:
                            with self.assertRaises(ValueError):
                                am.A_max(k, Ji, L1, L2, Jf)
        

    def test_A_max_equivalent_B_times_F_J_integral_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(a); integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                if int(L1) == int(L2):

                    try:
                        Ak_max = am.A_max(k, Ji, L1, L2, Jf)
                        if k == 2:
                            B2 = am.calc_B(k,Ji)
                            F2 = am.calc_F(k,Jf,L1,L2,Ji)
                            assert Ak_max == pytest.approx(B2*F2, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                        elif k == 4:
                            B4 = am.calc_B(k,Ji)
                            F4 = am.calc_F(k,Jf,L1,L2,Ji)
                            assert Ak_max == pytest.approx(B4*F4, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                    except ValueError:
                        if not B2 or F2 or B4 or F4:
                            with self.assertRaises(ValueError):
                                am.A_max(k, Ji, L1, L2, Jf)

    def test_A_max_equivalent_B_times_F_J_halfint_and_returns_float(self):
        # Compare calculated A_max with BF from Table 2(a); integral J
        # Ensure floating-point data type gets returned
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_orders = [2,4]
        for k in k_orders:
            for row in ldf:
                Ji = row[0]
                Jf = row[1]
                L1 = row[2]
                L2 = row[3]

                if int(L1) == int(L2):

                    try:
                        Ak_max = am.A_max(k, Ji, L1, L2, Jf)
                        if k == 2:
                            B2 = am.calc_B(k,Ji)
                            F2 = am.calc_F(k,Jf,L1,L2,Ji)
                            assert Ak_max == pytest.approx(B2*F2, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                        elif k == 4:
                            B4 = am.calc_B(k,Ji)
                            F4 = am.calc_F(k,Jf,L1,L2,Ji)
                            assert Ak_max == pytest.approx(B4*F4, abs=0.00001)
                            self.assertIsInstance(Ak_max, float)
                    except ValueError:
                        if not B2 or F2 or B4 or F4:
                            with self.assertRaises(ValueError):
                                am.A_max(k, Ji, L1, L2, Jf)

    # Type tests

    def test_A_max_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.A_max()
        # Too few args
        with self.assertRaises(TypeError):
            am.A_max(2,2,2,2)
        # Too many args
        with self.assertRaises(TypeError):
            # default 6th argument dg=0
            am.A_max(2,2,2,2,0,0.25,2)

    def test_A_max_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.A_max("2","2","2","2","0")
        with self.assertRaises(TypeError):
            am.A_max("2","2","2","2","0","0.25")
            
    def test_A_max_raises_TypeError_with_k_odd_arg(self):
        # A_max method tries to return B*F
        # Odd-k attempts NoneType * NoneType
        with self.assertRaises(TypeError):
            am.A_max(1,2,2,2,0)
        with self.assertRaises(TypeError):
            am.A_max(3,2,2,2,0)
        with self.assertRaises(TypeError):
            am.A_max(5,2,2,2,0)

