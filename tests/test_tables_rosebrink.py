import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

def fill_list_with_zeros(my_list, desired_length=7):
    current_length = len(my_list)
    if current_length < desired_length:
        # Calculate how many zeros are needed
        zeros_to_add = desired_length - current_length
        # Extend the list with the required number of zeros
        my_list.extend([0] * zeros_to_add)
    return my_list


class RoseBrinkTestsR(unittest.TestCase):

    __doc__="""Unit tests for table generation and manipulation methods 
    corresponding to the data presented in Rose and Brink's paper [2].

    R-coefficients

    References:
    [2] H.J. Rose, D.M. Brink, Rev. Mod. Phys. 39, 306 (1967).
    """

    # Table R tests
    # Integral J
    def test_get_tableRa_returns_DataFrame_of_length_61(self):
        df = am.get_tableRa()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 61
    
    def test_get_tableRa_check_data_types(self):
        df = am.get_tableRa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 61
        for v in ldf:
            assert len(v) == 16
            integral_Ji = int(2*v[0])
            assert integral_Ji % 2 == 0
            self.assertIsInstance(integral_Ji, int)

            integral_Jf = int(2*v[1])
            assert integral_Jf % 2 == 0
            self.assertIsInstance(integral_Jf, int)

            integral_L1 = int(2*v[2])
            assert integral_L1 % 2 == 0
            self.assertIsInstance(integral_L1, int)

            integral_L2 = int(2*v[3])
            assert integral_L2 % 2 == 0
            self.assertIsInstance(integral_L2, int)

            for i in range(4,16):
                self.assertIsInstance(v[i],float)

    def test_get_tableRa_check_return_against_calc_R(self):
        df = am.get_tableRa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            Ji = int(v[0])
            Jf = int(v[1])
            L1 = int(v[2])
            L2 = int(v[3])

            # k = 2
            try:
                assert v[4] == pytest.approx(am.calc_R(2,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[4]

            try:
                assert v[5] == pytest.approx(am.calc_R(2,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[5]

            try:
                assert v[6] == pytest.approx(am.calc_R(2,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[6]

            # k = 4
            try:
                assert v[7] == pytest.approx(am.calc_R(4,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[7]

            try:
                assert v[8] == pytest.approx(am.calc_R(4,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[8]

            try:
                assert v[9] == pytest.approx(am.calc_R(4,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[9]

            # k = 6
            try:
                assert v[10] == pytest.approx(am.calc_R(6,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[10]

            try:
                assert v[11] == pytest.approx(am.calc_R(6,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[11]

            try:
                assert v[12] == pytest.approx(am.calc_R(6,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[12]

            # k = 8
            try:
                assert v[13] == pytest.approx(am.calc_R(8,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[13]

            try:
                assert v[14] == pytest.approx(am.calc_R(8,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[14]

            try:
                assert v[15] == pytest.approx(am.calc_R(8,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[15]

    def test_get_tableRa_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tableRa(ARG)

    # Half-integral J
    def test_get_tableRb_returns_DataFrame_of_length_54(self):
        df = am.get_tableRb()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 54
    
    def test_get_tableRb_check_data_types(self):
        df = am.get_tableRb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 54
        for v in ldf:
            assert len(v) == 16
            half_integral_Ji = int(2*v[0])
            assert half_integral_Ji % 2 == 1
            self.assertIsInstance(v[0], float)

            half_integral_Jf = int(2*v[1])
            assert half_integral_Jf % 2 == 1
            self.assertIsInstance(v[1], float)

            integral_L1 = int(2*v[2])
            assert integral_L1 % 2 == 0
            self.assertIsInstance(integral_L1, int)

            integral_L2 = int(2*v[3])
            assert integral_L2 % 2 == 0
            self.assertIsInstance(integral_L2, int)

            for i in range(4,16):
                self.assertIsInstance(v[i],float)

    def test_get_tableRb_check_return_against_calc_R(self):
        df = am.get_tableRb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            Ji = float(v[0])
            Jf = float(v[1])
            L1 = int(v[2])
            L2 = int(v[3])

            # k = 2
            try:
                assert v[4] == pytest.approx(am.calc_R(2,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[4]

            try:
                assert v[5] == pytest.approx(am.calc_R(2,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[5]

            try:
                assert v[6] == pytest.approx(am.calc_R(2,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[6]

            # k = 4
            try:
                assert v[7] == pytest.approx(am.calc_R(4,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[7]

            try:
                assert v[8] == pytest.approx(am.calc_R(4,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[8]

            try:
                assert v[9] == pytest.approx(am.calc_R(4,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[9]

            # k = 6
            try:
                assert v[10] == pytest.approx(am.calc_R(6,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[10]

            try:
                assert v[11] == pytest.approx(am.calc_R(6,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[11]

            try:
                assert v[12] == pytest.approx(am.calc_R(6,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[12]

            # k = 8
            try:
                assert v[13] == pytest.approx(am.calc_R(8,L1,L1,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[13]

            try:
                assert v[14] == pytest.approx(am.calc_R(8,L1,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[14]

            try:
                assert v[15] == pytest.approx(am.calc_R(8,L2,L2,Ji,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[15]

    def test_get_tableRb_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tableRb(ARG)


class RoseBrinkTestsU(unittest.TestCase):

    __doc__="""Unit tests for table generation and manipulation methods 
    corresponding to the data presented in Rose and Brink's paper [2].

    U-coefficients

    References:
    [2] H.J. Rose, D.M. Brink, Rev. Mod. Phys. 39, 306 (1967).
    """

    def test_get_tableU_returns_DataFrame_of_length_75(self):
        df = am.get_tableU()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 75

    def test_get_tableU_check_data_types(self):
        df = am.get_tableU()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 75
        for v in ldf:
            assert len(v) == 12
            
            self.assertIsInstance(v[0], float)
            self.assertIsInstance(v[1], float)

            integral_L1 = int(2*v[2])
            assert integral_L1 % 2 == 0
            self.assertIsInstance(integral_L1, int)

            integral_L2 = int(2*v[3])
            assert integral_L2 % 2 == 0
            self.assertIsInstance(integral_L2, int)

            for i in range(4,12):
                self.assertIsInstance(v[i],float)

    def test_get_tableU_check_return_against_calc_u(self):
        df = am.get_tableU()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            Ji = float(v[0])
            Jf = float(v[1])
            L1 = int(v[2])
            L2 = int(v[3])

            # k = 2
            try:
                assert v[4] == pytest.approx(am.calc_u(2,Ji,L1,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[4]

            try:
                assert v[5] == pytest.approx(am.calc_u(2,Ji,L2,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[5]

            # k = 4
            try:
                assert v[6] == pytest.approx(am.calc_u(4,Ji,L1,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[6]

            try:
                assert v[7] == pytest.approx(am.calc_u(4,Ji,L2,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[7]

            # k = 6
            try:
                assert v[8] == pytest.approx(am.calc_u(6,Ji,L1,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[8]

            try:
                assert v[9] == pytest.approx(am.calc_u(6,Ji,L2,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[9]

            # k = 8
            try:
                assert v[10] == pytest.approx(am.calc_u(8,Ji,L1,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[10]

            try:
                assert v[11] == pytest.approx(am.calc_u(8,Ji,L2,Jf), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[11]

    def test_get_tableU_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tableU(ARG)


class RoseBrinkTestsU(unittest.TestCase):

    __doc__="""Unit tests for table generation and manipulation methods 
    corresponding to the data presented in Rose and Brink's paper [2].

    S-coefficients

    References:
    [2] H.J. Rose, D.M. Brink, Rev. Mod. Phys. 39, 306 (1967).
    """

    # Integral J and s
    def test_get_tableSa_returns_DataFrame_of_length_105(self):
        df = am.get_tableSa()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 105

    def test_get_tableSa_check_data_types(self):
        df = am.get_tableSa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 105
        for v in ldf:
            assert len(v) == 9
            
            self.assertIsInstance(int(v[0]), int)
            self.assertIsInstance(int(v[1]), int)
            self.assertIsInstance(float(v[2]), float)
            self.assertIsInstance(float(v[3]), float)

            for i in range(4,9):
                self.assertIsInstance(v[i],float)

    def test_get_tableSa_check_return_against_calc_S(self):
        df = am.get_tableSa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            l1 = float(v[0])
            l2 = float(v[1])
            J = float(v[2])
            s = float(v[3])

            # k = 0
            try:
                assert v[4] == pytest.approx(am.calc_S(l1,l2,J,s,0), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[4]

            # k = 2
            try:
                assert v[5] == pytest.approx(am.calc_S(l1,l2,J,s,2), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[5]

            # k = 4
            try:
                assert v[6] == pytest.approx(am.calc_S(l1,l2,J,s,4), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[6]

            # k = 6
            try:
                assert v[7] == pytest.approx(am.calc_S(l1,l2,J,s,6), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[7]

            # k = 8
            try:
                assert v[8] == pytest.approx(am.calc_S(l1,l2,J,s,8), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[8]

    def test_get_tableSa_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tableSa(ARG)

    # Half-integral J and s
    def test_get_tableSb_returns_DataFrame_of_length_105(self):
        df = am.get_tableSb()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 83

    def test_get_tableSb_check_data_types(self):
        df = am.get_tableSb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 83
        for v in ldf:
            assert len(v) == 9
            
            self.assertIsInstance(int(v[0]), int)
            self.assertIsInstance(int(v[1]), int)
            assert (2*float(v[2])) % 2 == 1
            self.assertIsInstance(float(v[2]), float)
            assert (2*float(v[3])) % 2 == 1
            self.assertIsInstance(float(v[3]), float)

            for i in range(4,9):
                self.assertIsInstance(v[i],float)

    def test_get_tableSb_check_return_against_calc_S(self):
        df = am.get_tableSb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            l1 = float(v[0])
            l2 = float(v[1])
            J = float(v[2])
            s = float(v[3])

            # k = 0
            try:
                assert v[4] == pytest.approx(am.calc_S(l1,l2,J,s,0), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[4]

            # k = 2
            try:
                assert v[5] == pytest.approx(am.calc_S(l1,l2,J,s,2), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[5]

            # k = 4
            try:
                assert v[6] == pytest.approx(am.calc_S(l1,l2,J,s,4), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[6]

            # k = 6
            try:
                assert v[7] == pytest.approx(am.calc_S(l1,l2,J,s,6), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[7]

            # k = 8
            try:
                assert v[8] == pytest.approx(am.calc_S(l1,l2,J,s,8), abs=0.0001)
            except:
                with self.assertRaises(TypeError):
                    assert not v[8]

    def test_get_tableSb_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tableSb(ARG)
            

class RoseBrinkTestsP(unittest.TestCase):

    __doc__="""Unit tests for table generation and manipulation methods 
    corresponding to the data presented in Rose and Brink's paper [2].

    `rho`-coefficients

    References:
    [2] H.J. Rose, D.M. Brink, Rev. Mod. Phys. 39, 306 (1967).
    """

    # Integral J
    def test_get_tablePa_returns_DataFrame_of_length_42(self):
        df = am.get_tablePa()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 42

    def test_get_tablePa_check_data_types(self):
        df = am.get_tablePa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 42
        for v in ldf:
            assert len(v) == 9

            self.assertIsInstance(int(v[0]), int)
            assert int(2*v[0]) % 2 == 0
            self.assertIsInstance(int(v[1]), int)
            assert int(2*v[1]) % 2 == 0

            for i in range(2,9):
                self.assertIsInstance(v[i],float)

    def test_get_tablePa_check_return_against_calc_p(self):
        df = am.get_tablePa()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            J = int(v[0])
            k = int(v[1])

            msubs = am.calc_p(k,J)

            if len(msubs) < 7:
                msubs = fill_list_with_zeros(msubs)

            assert v[2] == pytest.approx(msubs[0], abs=0.0001)
            assert v[3] == pytest.approx(msubs[1], abs=0.0001)
            assert v[4] == pytest.approx(msubs[2], abs=0.0001)
            assert v[5] == pytest.approx(msubs[3], abs=0.0001)
            assert v[6] == pytest.approx(msubs[4], abs=0.0001)
            assert v[7] == pytest.approx(msubs[5], abs=0.0001)
            assert v[8] == pytest.approx(msubs[6], abs=0.0001)

    def test_get_tablePa_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tablePa(ARG)

    # Half-integral J
    def test_get_tablePb_returns_DataFrame_of_length_38(self):
        df = am.get_tablePb()
        self.assertIsInstance(df, pd.core.frame.DataFrame)
        assert len(df) == 38

    def test_get_tablePa_check_data_types(self):
        df = am.get_tablePb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        assert len(ldf) == 38
        for v in ldf:
            assert len(v) == 9

            self.assertIsInstance(int(v[0]), int)
            assert int(2*v[0]) % 2 == 1
            self.assertIsInstance(int(v[1]), int)
            assert int(2*v[1]) % 2 == 0

            for i in range(2,9):
                self.assertIsInstance(v[i],float)

    def test_get_tablePa_check_return_against_calc_p(self):
        df = am.get_tablePb()
        ldf = df.values.tolist()
        self.assertIsInstance(ldf, list)
        for v in ldf:                
            J = float(v[0])
            k = int(v[1])

            msubs = am.calc_p(k,J)

            if len(msubs) < 7:
                msubs = fill_list_with_zeros(msubs)

            assert v[2] == pytest.approx(msubs[0], abs=0.0001)
            assert v[3] == pytest.approx(msubs[1], abs=0.0001)
            assert v[4] == pytest.approx(msubs[2], abs=0.0001)
            assert v[5] == pytest.approx(msubs[3], abs=0.0001)
            assert v[6] == pytest.approx(msubs[4], abs=0.0001)
            assert v[7] == pytest.approx(msubs[5], abs=0.0001)
            assert v[8] == pytest.approx(msubs[6], abs=0.0001)

    def test_get_tablePa_raises_TypeError_when_passing_args(self):
        ARG = 1
        with self.assertRaises(TypeError):
            am.get_tablePb(ARG)            
