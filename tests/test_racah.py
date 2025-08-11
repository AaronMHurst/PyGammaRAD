import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class RacahTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the Racah class for the 
    evaluation of Racha coefficients resulting from the coupling of a set of 
    angular momenta.

    Several tests have been written to confirm the results of evaluated 
    Racah coefficients used in the propagation of angular distribution 
    functions published in Tables 2(a) and 2(b) of Ref. [Yamazaki 1967].

    [Yamazaki 1967] T. Yamazaki, Nucl. Data Sheets A 3, 1 (1967).
    """

    # Tests confirming known evaluated Racah coefficients

    def test_racah_confirm_value1(self):
        W = am.racah(1,1,1,1,2,0)
        expected_W = 1/3
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value2(self):
        W = am.racah(1,1,1,1,2,1)
        expected_W = 1/6
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value3(self):
        W = am.racah(1,1,1,2,2,2)
        expected_W = 1/10
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value4(self):
        W = am.racah(5,5,2,2,2,4)
        expected_W = 0
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value5(self):
        W = am.racah(1.5,1.5,1,2,2,1.5)
        expected_W = np.sqrt(2)/10
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value6(self):
        W = am.racah(1.5,1.5,2,2,2,1.5)
        expected_W = 0
        assert W == pytest.approx(expected_W)



    def test_racah_confirm_value7(self):
        W = am.racah(1,1,1,1,2,2)
        expected_W = 1/30
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value8(self):
        W = am.racah(4,4,4,4,4,1)
        expected_W = -1/18
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value9(self):
        W = am.racah(1.5,1.5,1.5,1.5,2,2)
        expected_W = 3/20
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value10(self):
        W = am.racah(2.5,2.5,2.5,2.5,4,2)
        expected_W = 1/12
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value11(self):
        W = am.racah(3.5,3.5,3.5,3.5,4,2)
        expected_W = 1/24
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value12(self):
        W = am.racah(4.5,4.5,4.5,4.5,4,2)
        expected_W = -1/660
        assert W == pytest.approx(expected_W)

    def test_racah_confirm_value13(self):
        W = am.racah(15,15,18,18,2,3)
        expected_W = 0.02907
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value14(self):
        W = am.racah(15,15,18,18,4,3)
        expected_W = 0.02801
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value15(self):
        W = am.racah(15,15,2,3,2,17)
        expected_W = 0.04773
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value16(self):
        W = am.racah(15,15,2,3,4,17)
        expected_W = 0.01793
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value17(self):
        W = am.racah(15.5,15.5,18.5,18.5,2,3)
        expected_W = -0.02825
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value18(self):
        W = am.racah(15.5,15.5,18.5,18.5,4,3)
        expected_W = -0.02729
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value19(self):
        W = am.racah(14.5,14.5,2,3,2,12.5)
        expected_W = -0.04843
        assert W == pytest.approx(expected_W, abs=0.00001)

    def test_racah_confirm_value20(self):
        W = am.racah(14.5,14.5,2,3,4,12.5)
        expected_W = -0.02886
        assert W == pytest.approx(expected_W, abs=0.00001)        

    # Type tests
    
    def test_racah_returns_float_type(self):
        W = am.racah(1.5,1.5,1.5,1.5,2,2)
        self.assertIsInstance(W, float)

    def test_racah_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.racah()
        # Too few args
        with self.assertRaises(TypeError):
            am.racah(1.5,1.5,1.5)
        # Too many args
        with self.assertRaises(TypeError):
            am.racah(1.5,1.5,1.5,1.5,2,2,2.5)

    def test_racah_raises_TypeError_with_wrong_types_args(self):
        # Pass args as str types
        with self.assertRaises(TypeError):
            am.racah("1.5","1.5","1.5","1.5","2","2")
