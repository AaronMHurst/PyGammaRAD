import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

def fill_list_with_zeros(my_list, desired_length=7):
    """Method to fill list of length below a desired value with 
    [0.0000]-valued elements."""
    current_length = len(my_list)
    if current_length < desired_length:
        # Calculate how many zeros are needed
        zeros_to_add = desired_length - current_length
        # Extend the list with the required number of zeros
        my_list.extend([0.0000] * zeros_to_add)
    return my_list


class AngularDistributionPkTests(unittest.TestCase):

    __doc__="""This module contains unit tests comparing the angular 
    distribution function for the Sk-coefficient calculated in accordance with 
    Equation (3.63) from Rose and Brink's review article paper 

    [1967Ya05] to those obtained from 
    Tables 2a (integral J) and 2b (half-integral J) of the same reference 
    [1967Ya05].  The return types are also assessed.
    
    References:
    [1967Ya05] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    
    def test_calc_p_equal_to_tablePa_values_and_returns_float_integral_J(self):
        # Test to compare calculated `rho` coefficients with tabulated values.
        # Loop over integral spins from `rho`-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tablePa()
        ldf = df.values.tolist()

        for row in ldf:
            J = row[0]
            k = row[1]

            m0_tabulated = row[2]
            m1_tabulated = row[3]
            m2_tabulated = row[4]
            m3_tabulated = row[5]
            m4_tabulated = row[6]
            m5_tabulated = row[7]
            m6_tabulated = row[8]

            mk_table_list = [m0_tabulated, m1_tabulated, m2_tabulated,
                             m3_tabulated, m4_tabulated, m5_tabulated,
                             m6_tabulated]

            for i,m in enumerate(range(0,7)):
                
                try:
                    m_calculated = am.calc_p(k, J, m)
                    # Using "m<i> == row[i+2]" relationship
                    assert m_calculated == pytest.approx(row[i+2], abs=0.0001)
                    self.assertIsInstance(m_calculated, float)

                except TypeError:
                    if int(row[i+2]) == 0: 
                        with self.assertRaises(TypeError):
                            am.calc_p(k, J, m)


            # Test calling method for all m-substates
            mk_calc_list = am.calc_p(k, J)
            if len(mk_calc_list) < 7:
                mk_calc_list = fill_list_with_zeros(mk_calc_list)
                
            for i,mk_t in enumerate(mk_table_list):
                for j,mk_c in enumerate(mk_calc_list):
                    if i == j:
                        assert mk_t == pytest.approx(mk_c, abs=0.0001)
                        self.assertIsInstance(mk_t, float)
                        self.assertIsInstance(mk_c, float)
                        
    def test_calc_p_equal_to_tablePb_values_and_returns_float_halfint_J(self):
        # Test to compare calculated `rho` coefficients with tabulated values.
        # Loop over half-integral spins from `rho`-coefficient Table 
        # Test to ensure return of floating-point value type
        df = am.get_tablePb()
        ldf = df.values.tolist()

        for row in ldf:
            J = row[0]
            k = row[1]

            m1_2_tabulated = row[2]
            m3_2_tabulated = row[3]
            m5_2_tabulated = row[4]
            m7_2_tabulated = row[5]
            m9_2_tabulated = row[6]
            m11_2_tabulated = row[7]
            m13_2_tabulated = row[8]

            mk_table_list = [m1_2_tabulated, m3_2_tabulated, m5_2_tabulated,
                             m7_2_tabulated, m9_2_tabulated, m11_2_tabulated,
                             m13_2_tabulated]

            for i,m in enumerate(range(1,14,2)):
                if m % 2 == 1:
                    m = m/2.0
                    try:
                        m_calculated = am.calc_p(k, J, m)
                        # Using "m<i> == row[i+2]" relationship
                        assert m_calculated == pytest.approx(row[i+2], abs=0.0001)
                        self.assertIsInstance(m_calculated, float)

                    except TypeError:
                        if int(row[i+2]) == 0: 
                            with self.assertRaises(TypeError):
                                am.calc_p(k, J, m)


            # Test calling method for all m-substates
            mk_calc_list = am.calc_p(k, J)
            if len(mk_calc_list) < 7:
                mk_calc_list = fill_list_with_zeros(mk_calc_list)
                
            for i,mk_t in enumerate(mk_table_list):
                for j,mk_c in enumerate(mk_calc_list):
                    if i == j:
                        assert mk_t == pytest.approx(mk_c, abs=0.0001)
                        self.assertIsInstance(mk_t, float)
                        self.assertIsInstance(mk_c, float)
                        

    # Type tests
    
    def test_calc_p_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.calc_p()
        # Too few args
        with self.assertRaises(TypeError):
            am.calc_p(2)

    def test_calc_p_returns_NoneType_with_too_many_args(self):
        # Too many args
        p = am.calc_p(2,1.5,0.5,3)
        self.assertIsNone(p)
        p = am.calc_p(8,10,5,7)
        self.assertIsNone(p)
        p = am.calc_p(6,15,0,4,2,0.5)
        self.assertIsNone(p)

            
    def test_calc_p_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.calc_p("8","10")
        with self.assertRaises(TypeError):
            am.calc_p("2","1.5","0.5")

    def test_calc_p_returns_NoneType_with_k_odd_arg(self):
        # Pass odd-k args
        p = am.calc_p(1,10)
        self.assertIsNone(p)
        p = am.calc_p(3,10,2)
        self.assertIsNone(p)
        p = am.calc_p(5,2.5,1.5)
        self.assertIsNone(p)
        p = am.calc_p(7,7.5)
        self.assertIsNone(p)
        p = am.calc_p(9,9.5,7.5)
        self.assertIsNone(p)
