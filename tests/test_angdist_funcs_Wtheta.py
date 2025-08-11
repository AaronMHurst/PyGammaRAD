import pytest
import unittest
import numpy as np
import pandas as pd

import PyGammaRAD as pg
am = pg.AngularMomentum()

class AngularDistributionWthetaTests(unittest.TestCase):

    __doc__="""This module contains unit tests to help validate the `dist_W` 
    method based on Eq. (2) in Yamazaki's paper [1].

    References:
    [1] T. Yamazaki Nucl. Data Sect. A, Vol. 3, Num. 1 (1967).
    """
    def test_dist_W_returns_tuple_containing_numpy_arrays_of_floats(self):
        # Default range
        W20 = am.dist_W(am.A_max(2,2,2,2,0),am.A_max(4,2,2,2,0))
        self.assertIsInstance(W20, tuple)
        assert len(W20) == 2
        assert len(W20[0]) == 1800
        for w in W20[0]:
            self.assertIsInstance(w, float)
        assert len(W20[1]) == 1800
        for w in W20[1]:
            self.assertIsInstance(w, float)

        user_range = np.linspace(0,90,90)
        assert len(user_range) == 90
        W20_user = am.dist_W(am.A_max(2,2,2,2,0),am.A_max(4,2,2,2,0),user_range)
        self.assertIsInstance(W20_user, tuple)
        assert len(W20_user) == 2
        assert len(W20_user[0]) == len(user_range)
        for w_user in W20_user[0]:
            self.assertIsInstance(w_user, float)
        assert len(W20_user[1]) == len(user_range)
        for w_user in W20_user[1]:
            self.assertIsInstance(w_user, float)
           
    # Type tests

    def test_dist_W_raises_TypeError_with_wrong_number_args(self):
        # No args, i.e., too few since only 1 arg required
        with self.assertRaises(TypeError):
            am.dist_W()

    def test_dist_W_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.dist_W("0.5","-0.4")
