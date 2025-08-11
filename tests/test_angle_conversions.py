import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngleConversionTests(unittest.TestCase):

    __doc__="""This module contains unit tests to help validate the conversion 
    methods of radians-to-degrees and vice versa contained within the 
    `Legendre` class.
    """

    ##### Degrees-To-Radians: `deg2rad` #####
    
    def test_deg2rad_returns_array_of_floats_from_array_or_list(self):
        # default
        deg = am.deg2rad()
        assert len(deg) == 1800
        self.assertIsInstance(deg, np.ndarray)
        for d in deg:
            self.assertIsInstance(d, float)

        # user-defined list
        theta_list = [31,12,43,4,4,74,24,4,76,5,2,35,3,2,78,5,10,85]
        self.assertIsInstance(theta_list, list)
        assert len(theta_list) == 18
        deg2 = am.deg2rad(theta_list)
        assert len(deg2) == len(theta_list)
        self.assertIsInstance(deg2, np.ndarray)
        for d2 in deg2:
            self.assertIsInstance(d2, float)

        # user-defined array
        theta_array = np.array([31,12,43,4,4,74,24,4,76,5,2,35,3,2,78,5,10,85])
        self.assertIsInstance(theta_array, np.ndarray)
        assert len(theta_array) == 18
        deg3 = am.deg2rad(theta_array)
        assert len(deg3) == len(theta_array)
        self.assertIsInstance(deg3, np.ndarray)
        for d3 in deg3:
            self.assertIsInstance(d3, float)

    def test_deg2rad_returns_single_value_float_from_float_or_integer(self):
        # single-value integer
        i = 43
        self.assertIsInstance(i, int)
        deg_i = am.deg2rad(i)
        self.assertIsInstance(deg_i, float)

        # single-value float
        f = 19.43
        self.assertIsInstance(f, float)
        deg_f = am.deg2rad(f)
        self.assertIsInstance(deg_f, float)

    def test_deg2rad_conversion_value_check(self):
        # Check conversion: theta_rad = (theta_deg/180) * pi
        theta_deg = [0,45,90,135,180]
        assert am.deg2rad(theta_deg[0]) == pytest.approx(0)
        assert am.deg2rad(theta_deg[1]) == pytest.approx(np.pi/4)
        assert am.deg2rad(theta_deg[2]) == pytest.approx(np.pi/2)
        assert am.deg2rad(theta_deg[3]) == pytest.approx(3*np.pi/4)
        assert am.deg2rad(theta_deg[4]) == pytest.approx(np.pi)
        
    # Type tests

    def test_deg2rad_raises_TypeError_with_wrong_args(self):
        theta_deg = [0,45,90,135,180]
        # Too many args
        with self.assertRaises(TypeError):
            am.deg2rad(45,90)
        with self.assertRaises(TypeError):
            am.deg2rad(theta_deg,45)
        # Wrong args types
        with self.assertRaises(TypeError):
            am.deg2rad("theta_deg")
        with self.assertRaises(TypeError):
            am.deg2rad("45")

    ##### Radians-To-Degrees: `rad2deg` #####
    
    def test_rad2deg_returns_array_of_floats_from_array_or_list(self):
        # default
        rad = am.rad2deg()
        assert len(rad) == 1800
        self.assertIsInstance(rad, np.ndarray)
        for r in rad:
            self.assertIsInstance(r, float)

        
        # user-defined list
        theta_list = [3.1,1.2,0.43,0.4,0.4,0.74,2.4,0.4,0.76,0.5,2.0,0.35,3.0,2.0,0.78,0.5,1.0,0.85]
        self.assertIsInstance(theta_list, list)
        assert len(theta_list) == 18
        rad2 = am.rad2deg(theta_list)
        assert len(rad2) == len(theta_list)
        self.assertIsInstance(rad2, np.ndarray)
        for r2 in rad2:
            self.assertIsInstance(r2, float)

        
        # user-defined array
        theta_array = np.array([3.1,1.2,0.43,0.4,0.4,0.74,2.4,0.4,0.76,0.5,2.0,0.35,3.0,2.0,0.78,0.5,1.0,0.85])
        self.assertIsInstance(theta_array, np.ndarray)
        assert len(theta_array) == 18
        rad3 = am.rad2deg(theta_array)
        assert len(rad3) == len(theta_array)
        self.assertIsInstance(rad3, np.ndarray)
        for r3 in rad3:
            self.assertIsInstance(r3, float)

    def test_rad2deg_returns_single_value_float_from_float_or_integer(self):
        # single-value integer
        i = 1
        self.assertIsInstance(i, int)
        rad_i = am.rad2deg(i)
        self.assertIsInstance(rad_i, float)

        # single-value float
        f = 0.1943
        self.assertIsInstance(f, float)
        rad_f = am.rad2deg(f)
        self.assertIsInstance(rad_f, float)

    def test_rad2deg_conversion_value_check(self):
        # Check conversion: theta_rad = (theta_rad/pi) * 180
        theta_rad = [0.0,np.pi/4.,np.pi/2.,3*np.pi/4.,np.pi]
        assert int(am.rad2deg(theta_rad[0])) == 0
        assert int(am.rad2deg(theta_rad[1])) == 45
        assert int(am.rad2deg(theta_rad[2])) == 90
        assert int(am.rad2deg(theta_rad[3])) == 135
        assert int(am.rad2deg(theta_rad[4])) == 180
        
    # Type tests

    def test_rad2deg_raises_TypeError_with_wrong_args(self):
        theta_rad = [0.0,np.pi/4.,np.pi/2.,3.*np.pi/4.,np.pi]
        # Too many args
        with self.assertRaises(TypeError):
            am.rad2deg(np.pi/4.,np.pi/2.)
        with self.assertRaises(TypeError):
            am.rad2deg(theta_rad,np.pi/4)
        # Wrong args types
        with self.assertRaises(TypeError):
            am.rad2deg("theta_rad")
        with self.assertRaises(TypeError):
            am.rad2deg("0.8787891")            
    
