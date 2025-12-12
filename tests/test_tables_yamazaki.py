import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class YamazakiTableTests(unittest.TestCase):

    __doc__="""Unit tests for table generation and manipulation methods 
    corresponding to the data presented in Yamazaki's paper [1].

    References:
    [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """

    # Table 1 tests
    
    def test_get_table1_returns_DataFrame_of_length_40(self):
        df = am.get_table1()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 40

    def test_get_table1_ensure_integral_spins(self):
        df = am.get_table1(0)
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 20
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 20
        # Check for integral spins; m=0 argument
        for i,v in enumerate(ldf):
            # Check each row has 5 elements
            assert len(v) == 5
            integral_J = int(2*v[0])
            self.assertIsInstance(integral_J, int)
            assert integral_J%2 == 0
            integral_m = int(2*v[1])
            self.assertIsInstance(integral_m, int)
            assert integral_m == 0
            assert integral_m%2 == 0

            # Check remaining datatypes
            if i==2 or i==3 or i==4:
                self.assertIsInstance(v[i], float)

    def test_get_table1_ensure_half_integral_spins(self):
        df = am.get_table1(0.5)
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 20
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 20
        # Check for half-integral spins; m=0.5 argument
        for i,v in enumerate(ldf):
            # Check each row has 5 elements
            assert len(v) == 5
            half_integral_J = int(2*v[0])
            self.assertIsInstance(half_integral_J, int)
            assert half_integral_J%2 == 1
            half_integral_m = int(2*v[1])
            self.assertIsInstance(half_integral_m, int)
            assert half_integral_m == 1
            assert half_integral_m%2 == 1

            # Check remaining datatypes
            if i==2 or i==3 or i==4:
                self.assertIsInstance(v[i], float)

    def test_get_B_J_equals_8_values_from_table1_check(self):
        k_order = [2,4,6]
        J = 8
        for k in k_order:
            Bk = am.get_B(k,J)
            self.assertIsInstance(Bk, float)
            # Spot-check Bk values for J=8
            if k==2:
                assert Bk == pytest.approx(-1.12390, abs=0.00001)
            elif k==4:
                assert Bk == pytest.approx(1.14531, abs=0.00001)
            else:
                assert Bk == pytest.approx(-1.17175, abs=0.00001)


    def test_get_B_always_returns_float(self):
        k_order = [2,4,6]
        for J in range(1,21):
            # Integral J
            for k in k_order:
                Bk = am.get_B(k,J)
                self.assertIsInstance(Bk, float)
            # Half-integral J
            J = J+0.5
            for k in k_order:
                Bk = am.get_B(k,J)
                self.assertIsInstance(Bk, float)
    
    def test_get_B_returns_NoneType_when_k_or_J_out_of_range(self):
        # Test k values outside of range
        Bk = am.get_B(1,5)
        self.assertIsNone(Bk)
        Bk = am.get_B(3,5)
        self.assertIsNone(Bk)
        Bk = am.get_B(7,5)
        self.assertIsNone(Bk)
        Bk = am.get_B(0,2)
        self.assertIsNone(Bk)

        k_order = [2,4,6]
        # Test integral J values outside of range
        for J in range(21,30):
            for k in k_order:
                Bk = am.get_B(k,J)
                self.assertIsNone(Bk)

        # Test half-integral J values outside of range
        # J = 1/2
        for k in k_order:
            Bk = am.get_B(k,0.5)
            self.assertIsNone(Bk)

        # J > 41/2
        for J in range(21,30):
            J = J+0.5
            for k in k_order:
                Bk = am.get_B(k,J)
                self.assertIsNone(Bk)
        
    def test_get_B_raises_TypeError_with_wrong_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_B()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_B(2)
        # Too many args
        with self.assertRaises(TypeError):
            am.calc_B(2,5,10)
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_B("2","5")

    def test_get_row_table1_returns_list_of_3_floats(self):
        for J in range(1,21):
            # Integral spin
            row_tab1 = am.get_row_table1(J)
            self.assertIsInstance(row_tab1, list)
            assert len(row_tab1) == 3
            for c in row_tab1:
                self.assertIsInstance(c, float)

            # Half-inetgral spin
            J = J + 0.5
            row_tab1 = am.get_row_table1(J)
            self.assertIsInstance(row_tab1, list)
            assert len(row_tab1) == 3
            for c in row_tab1:
                self.assertIsInstance(c, float)

    def test_get_row_table1_spot_check_data(self):
        J = 11/2
        row_tab1 = am.get_row_table1(J)
        row_tab1[0] == pytest.approx(-1.11144, abs=0.00001)
        row_tab1[1] == pytest.approx(1.10240, abs=0.00001)
        row_tab1[2] == pytest.approx(-1.07749, abs=0.00001)

        J = 19
        row_tab1 = am.get_row_table1(J)
        row_tab1[0] == pytest.approx(-1.11914, abs=0.00001)
        row_tab1[1] == pytest.approx(1.12873, abs=0.00001)
        row_tab1[2] == pytest.approx(-1.13465, abs=0.00001)
                
    # Table 2(a) tests

    def test_get_table2a_returns_DataFrame_of_length_246(self):
        df = am.get_table2a()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 246
    
    def test_get_table2a_contains_integral_spins(self):
        df = am.get_table2a()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 246
        # Check Ji and Jf values are integral
        for i,v in enumerate(ldf):
            # Check each row has 10 elements
            assert len(v) == 10
            Ji_integral = int(2*v[0])
            Jf_integral = int(2*v[1])
            self.assertIsInstance(Ji_integral, int)
            self.assertIsInstance(Jf_integral, int)
            assert Ji_integral%2 == 0
            assert Jf_integral%2 == 0

            # Check remaining datatypes
            if i==2 or i==3:
                L = int(v[i])
                self.assertIsInstance(L, int)
            elif i>=4 and i<10:
                self.assertIsInstance(v[i], float)
            
    def test_get_row_table2_a_check_data(self):
        df = am.get_table2a()
        ldf = df.values.tolist()
        k_order = [2,4]
        for row in ldf:
            Ji = row[0]
            Jf = row[1]
            L1 = int(row[2])
            L2 = int(row[3])

            coeffs = am.get_row_table2(Ji,Jf,L1,L2)
            self.assertIsInstance(coeffs, list)
            assert len(coeffs) == 6
            for c in coeffs:
                self.assertIsInstance(c, float)

            coeffs_F = am.get_row_table2(Ji,Jf,L1,L2,coeff='F')
            self.assertIsInstance(coeffs_F, list)
            assert len(coeffs_F) == 2
            for cF in coeffs_F:
                self.assertIsInstance(cF, float)
                
            coeffs_BF = am.get_row_table2(Ji,Jf,L1,L2,coeff='BF')
            self.assertIsInstance(coeffs_BF, list)
            assert len(coeffs_BF) == 2
            for cBF in coeffs_BF:
                self.assertIsInstance(cBF, float)
                
            coeffs_U = am.get_row_table2(Ji,Jf,L1,L2,coeff='U')
            self.assertIsInstance(coeffs_U, list)
            assert len(coeffs_U) == 2
            for cU in coeffs_U:
                self.assertIsInstance(cU, float)

            for k in k_order:
                coeffs_k = am.get_row_table2(Ji,Jf,L1,L2,k)
                self.assertIsInstance(coeffs_k, list)
                assert len(coeffs_k) == 3
                for c_k in coeffs_k:
                    self.assertIsInstance(c_k, float)

                coeff_Fk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='F')
                self.assertIsInstance(coeff_Fk, float)
                coeff_BkFk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='BF')
                self.assertIsInstance(coeff_BkFk, float)
                coeff_Uk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='U')
                self.assertIsInstance(coeff_Uk, float)
                
    # Table 2(b) tests

    def test_get_table2b_returns_DataFrame_of_length_250(self):
        df = am.get_table2b()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 250        

    def test_get_table2b_contains_half_integral_spins(self):
        df = am.get_table2b()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 250
        # Check Ji and Jf values are half-integral
        for i,v in enumerate(ldf):
            # Check each row has 10 elements
            assert len(v) == 10
            Ji_half_integral = int(2*v[0])
            Jf_half_integral = int(2*v[1])
            self.assertIsInstance(Ji_half_integral, int)
            self.assertIsInstance(Jf_half_integral, int)
            assert Ji_half_integral%2 == 1
            assert Jf_half_integral%2 == 1

            # Check remaining datatypes
            if i==2 or i==3:
                L = int(v[i])
                self.assertIsInstance(L, int)
            elif i>=4 and i<10:
                self.assertIsInstance(v[i], float)

    def test_get_row_table2_b_check_data(self):
        df = am.get_table2b()
        ldf = df.values.tolist()
        k_order = [2,4]
        for row in ldf:
            Ji = row[0]
            Jf = row[1]
            L1 = int(row[2])
            L2 = int(row[3])

            coeffs = am.get_row_table2(Ji,Jf,L1,L2)
            self.assertIsInstance(coeffs, list)
            assert len(coeffs) == 6
            for c in coeffs:
                self.assertIsInstance(c, float)

            coeffs_F = am.get_row_table2(Ji,Jf,L1,L2,coeff='F')
            self.assertIsInstance(coeffs_F, list)
            assert len(coeffs_F) == 2
            for cF in coeffs_F:
                self.assertIsInstance(cF, float)
                
            coeffs_BF = am.get_row_table2(Ji,Jf,L1,L2,coeff='BF')
            self.assertIsInstance(coeffs_BF, list)
            assert len(coeffs_BF) == 2
            for cBF in coeffs_BF:
                self.assertIsInstance(cBF, float)
                
            coeffs_U = am.get_row_table2(Ji,Jf,L1,L2,coeff='U')
            self.assertIsInstance(coeffs_U, list)
            assert len(coeffs_U) == 2
            for cU in coeffs_U:
                self.assertIsInstance(cU, float)

            for k in k_order:
                coeffs_k = am.get_row_table2(Ji,Jf,L1,L2,k)
                self.assertIsInstance(coeffs_k, list)
                assert len(coeffs_k) == 3
                for c_k in coeffs_k:
                    self.assertIsInstance(c_k, float)

                coeff_Fk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='F')
                self.assertIsInstance(coeff_Fk, float)
                coeff_BkFk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='BF')
                self.assertIsInstance(coeff_BkFk, float)
                coeff_Uk = am.get_row_table2(Ji,Jf,L1,L2,k,coeff='U')
                self.assertIsInstance(coeff_Uk, float)

    # Type tests for Table 2 using wrong args

    def test_get_row_table2_returns_NoneType_out_of_range(self):
        # Outside J range
        r = am.get_row_table2(15.5,19.5,3,3)
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,2)
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,coeff="F")
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,4,coeff="F")
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,coeff="BF")
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,2,coeff="BF")
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,coeff="U")
        self.assertIsNone(r)
        r = am.get_row_table2(15.5,19.5,3,3,4,coeff="U")
        self.assertIsNone(r)

        # Outside k range
        r = am.get_row_table2(15.5,18.5,3,3,8)
        self.assertIsNone(r)

        # Outside L range
        r = am.get_row_table2(15.5,18.5,3,5)
        self.assertIsNone(r)

    def test_get_row_table2_raises_TypeError_with_wrong_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.get_row_table2()
        # Too few args
        with self.assertRaises(TypeError):
            am.get_row_table2(15.5,18.5,3)
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.get_row_table2("15.5","18.5","3","3","2")

