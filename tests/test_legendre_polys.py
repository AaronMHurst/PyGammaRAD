import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class LegendrePolynomialTests(unittest.TestCase):

    __doc__="""This module contains unit tests to help validate the Legendre 
    polnomial methods of the `Legendre` class.
    """

    def test_lpolys_all_return_default_array_of_floats(self):
        p0 = am.lpoly0()
        assert len(p0) == 1800
        self.assertIsInstance(p0, np.ndarray)
        for p in p0:
            self.assertIsInstance(p, float)

        p1 = am.lpoly1()
        assert len(p1) == 1800
        self.assertIsInstance(p1, np.ndarray)
        for p in p1:
            self.assertIsInstance(p, float)

        p2 = am.lpoly2()
        assert len(p2) == 1800
        self.assertIsInstance(p2, np.ndarray)
        for p in p2:
            self.assertIsInstance(p, float)

        p3 = am.lpoly3()
        assert len(p3) == 1800
        self.assertIsInstance(p3, np.ndarray)
        for p in p3:
            self.assertIsInstance(p, float)

        p4 = am.lpoly4()
        assert len(p4) == 1800
        self.assertIsInstance(p4, np.ndarray)
        for p in p4:
            self.assertIsInstance(p, float)

        p5 = am.lpoly5()
        assert len(p5) == 1800
        self.assertIsInstance(p5, np.ndarray)
        for p in p5:
            self.assertIsInstance(p, float)

        p6 = am.lpoly6()
        assert len(p6) == 1800
        self.assertIsInstance(p6, np.ndarray)
        for p in p6:
            self.assertIsInstance(p, float)

        p7 = am.lpoly7()
        assert len(p7) == 1800
        self.assertIsInstance(p7, np.ndarray)
        for p in p7:
            self.assertIsInstance(p, float)

        p8 = am.lpoly8()
        assert len(p8) == 1800
        self.assertIsInstance(p8, np.ndarray)
        for p in p8:
            self.assertIsInstance(p, float)

        p9 = am.lpoly9()
        assert len(p9) == 1800
        self.assertIsInstance(p9, np.ndarray)
        for p in p9:
            self.assertIsInstance(p, float)

        p10 = am.lpoly10()
        assert len(p10) == 1800
        self.assertIsInstance(p10, np.ndarray)
        for p in p10:
            self.assertIsInstance(p, float)

    def test_lpolys_all_return_user_range_array_of_floats(self):
        theta_range = np.linspace(45,135,90)
        assert len(theta_range) == 90

        P0 = am.lpoly0(theta_range)
        assert len(P0) == len(theta_range)
        self.assertIsInstance(P0, np.ndarray)
        for P in P0:
            self.assertIsInstance(P, float)

        P1 = am.lpoly1(theta_range)
        assert len(P1) == len(theta_range)
        self.assertIsInstance(P1, np.ndarray)
        for P in P1:
            self.assertIsInstance(P, float)

        P2 = am.lpoly2(theta_range)
        assert len(P2) == len(theta_range)
        self.assertIsInstance(P2, np.ndarray)
        for P in P2:
            self.assertIsInstance(P, float)

        P3 = am.lpoly3(theta_range)
        assert len(P3) == len(theta_range)
        self.assertIsInstance(P3, np.ndarray)
        for P in P3:
            self.assertIsInstance(P, float)

        P4 = am.lpoly4(theta_range)
        assert len(P4) == len(theta_range)
        self.assertIsInstance(P4, np.ndarray)
        for P in P4:
            self.assertIsInstance(P, float)

        P5 = am.lpoly5(theta_range)
        assert len(P5) == len(theta_range)
        self.assertIsInstance(P5, np.ndarray)
        for P in P5:
            self.assertIsInstance(P, float)

        P6 = am.lpoly6(theta_range)
        assert len(P6) == len(theta_range)
        self.assertIsInstance(P6, np.ndarray)
        for P in P6:
            self.assertIsInstance(P, float)

        P7 = am.lpoly7(theta_range)
        assert len(P7) == len(theta_range)
        self.assertIsInstance(P7, np.ndarray)
        for P in P7:
            self.assertIsInstance(P, float)

        P8 = am.lpoly8(theta_range)
        assert len(P8) == len(theta_range)
        self.assertIsInstance(P8, np.ndarray)
        for P in P8:
            self.assertIsInstance(P, float)

        P9 = am.lpoly9(theta_range)
        assert len(P9) == len(theta_range)
        self.assertIsInstance(P9, np.ndarray)
        for P in P9:
            self.assertIsInstance(P, float)

        P10 = am.lpoly10(theta_range)
        assert len(P10) == len(theta_range)
        self.assertIsInstance(P10, np.ndarray)
        for P in P10:
            self.assertIsInstance(P, float)
                    
    def test_lpolys_all_return_floats_when_passing_single_value(self):
        theta = 52.35
        self.assertIsInstance(theta, float)

        P0 = am.lpoly0(theta)
        self.assertIsInstance(P0, float)

        P1 = am.lpoly1(theta)
        self.assertIsInstance(P1, float)

        P2 = am.lpoly2(theta)
        self.assertIsInstance(P2, float)

        P3 = am.lpoly3(theta)
        self.assertIsInstance(P3, float)

        P4 = am.lpoly4(theta)
        self.assertIsInstance(P4, float)

        P5 = am.lpoly5(theta)
        self.assertIsInstance(P5, float)

        P6 = am.lpoly6(theta)
        self.assertIsInstance(P6, float)

        P7 = am.lpoly7(theta)
        self.assertIsInstance(P7, float)

        P8 = am.lpoly8(theta)
        self.assertIsInstance(P8, float)

        P9 = am.lpoly9(theta)
        self.assertIsInstance(P9, float)

        P10 = am.lpoly10(theta)
        self.assertIsInstance(P10, float)

    def test_lpolys_all_return_user_range_array_of_floats_from_list(self):
        theta_list = [31, 12, 43, 4, 4, 74, 24, 4, 76, 5, 2, 35, 3, 2, 78]
        assert len(theta_list) == 15
        self.assertIsInstance(theta_list, list)

        P0 = am.lpoly0(theta_list)
        assert len(P0) == len(theta_list)
        self.assertIsInstance(P0, np.ndarray)
        for P in P0:
            self.assertIsInstance(P, float)

        P1 = am.lpoly1(theta_list)
        assert len(P1) == len(theta_list)
        self.assertIsInstance(P1, np.ndarray)
        for P in P1:
            self.assertIsInstance(P, float)

        P2 = am.lpoly2(theta_list)
        assert len(P2) == len(theta_list)
        self.assertIsInstance(P2, np.ndarray)
        for P in P2:
            self.assertIsInstance(P, float)

        P3 = am.lpoly3(theta_list)
        assert len(P3) == len(theta_list)
        self.assertIsInstance(P3, np.ndarray)
        for P in P3:
            self.assertIsInstance(P, float)

        P4 = am.lpoly4(theta_list)
        assert len(P4) == len(theta_list)
        self.assertIsInstance(P4, np.ndarray)
        for P in P4:
            self.assertIsInstance(P, float)

        P5 = am.lpoly5(theta_list)
        assert len(P5) == len(theta_list)
        self.assertIsInstance(P5, np.ndarray)
        for P in P5:
            self.assertIsInstance(P, float)

        P6 = am.lpoly6(theta_list)
        assert len(P6) == len(theta_list)
        self.assertIsInstance(P6, np.ndarray)
        for P in P6:
            self.assertIsInstance(P, float)

        P7 = am.lpoly7(theta_list)
        assert len(P7) == len(theta_list)
        self.assertIsInstance(P7, np.ndarray)
        for P in P7:
            self.assertIsInstance(P, float)

        P8 = am.lpoly8(theta_list)
        assert len(P8) == len(theta_list)
        self.assertIsInstance(P8, np.ndarray)
        for P in P8:
            self.assertIsInstance(P, float)

        P9 = am.lpoly9(theta_list)
        assert len(P9) == len(theta_list)
        self.assertIsInstance(P9, np.ndarray)
        for P in P9:
            self.assertIsInstance(P, float)

        P10 = am.lpoly10(theta_list)
        assert len(P10) == len(theta_list)
        self.assertIsInstance(P10, np.ndarray)
        for P in P10:
            self.assertIsInstance(P, float)

    def test_lpolys_all_raise_TypeError_with_wrong_args(self):
        theta = np.array([10, 5, 85])
        assert len(theta) == 3
        self.assertIsInstance(theta, np.ndarray)

        with self.assertRaises(TypeError):
            am.lpoly0(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly0("theta")
        with self.assertRaises(TypeError):
            am.lpoly0("500.0")

        with self.assertRaises(TypeError):
            am.lpoly1(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly1("theta")
        with self.assertRaises(TypeError):
            am.lpoly1("500.0")

        with self.assertRaises(TypeError):
            am.lpoly2(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly2("theta")
        with self.assertRaises(TypeError):
            am.lpoly2("500.0")

        with self.assertRaises(TypeError):
            am.lpoly3(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly3("theta")
        with self.assertRaises(TypeError):
            am.lpoly3("500.0")

        with self.assertRaises(TypeError):
            am.lpoly4(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly4("theta")
        with self.assertRaises(TypeError):
            am.lpoly4("500.0")

        with self.assertRaises(TypeError):
            am.lpoly5(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly5("theta")
        with self.assertRaises(TypeError):
            am.lpoly5("500.0")

        with self.assertRaises(TypeError):
            am.lpoly6(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly6("theta")
        with self.assertRaises(TypeError):
            am.lpoly6("500.0")

        with self.assertRaises(TypeError):
            am.lpoly7(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly7("theta")
        with self.assertRaises(TypeError):
            am.lpoly7("500.0")

        with self.assertRaises(TypeError):
            am.lpoly8(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly8("theta")
        with self.assertRaises(TypeError):
            am.lpoly8("500.0")

        with self.assertRaises(TypeError):
            am.lpoly9(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly9("theta")
        with self.assertRaises(TypeError):
            am.lpoly9("500.0")

        with self.assertRaises(TypeError):
            am.lpoly10(theta,500)
        with self.assertRaises(TypeError):
            am.lpoly10("theta")
        with self.assertRaises(TypeError):
            am.lpoly10("500.0")

            
