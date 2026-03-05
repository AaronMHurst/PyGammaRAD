import pytest
import unittest
import numpy as np

import PyGammaRAD as pg
am = pg.AngularMomentum()

class GauntTests(unittest.TestCase):

    __doc__="""Unit tests for the `gaunt` method which invokes the associated 
    methods belonging to the `Wigner3j` class for the evaluation of Wigner-3j 
    symbols resulting from the coupling of a set of angular momenta.

    Several tests have been developed to confirm the numerical results of  
    published Gaunt coefficients tabulated in Table 1 of Ref. [Yukcu 2025].

    [Yukcu 2025] S.A. Yukcu, Can. J. Phys. 103, 321 (2025).
    """

    # Tests confirming tabulated Gaunt coefficients

    def test_gaunt_confirm_value1(self):
        G = am.gaunt(2,58,60,-1,-2,3)
        expected_G = -3*np.sqrt(16347695)/(34034*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value2(self):
        G = am.gaunt(6,15,15,0,0,0)
        expected_G = 1768*np.sqrt(13)/(35409*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value3(self):
        G = am.gaunt(10,10,12,9,3,-12)
        expected_G = -98*np.sqrt(6279)/(62031*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-15)

    def test_gaunt_confirm_value4(self):
        G = am.gaunt(12,15,5,2,3,-5)
        expected_G = 91*np.sqrt(36890)/(124062*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value5(self):
        G = am.gaunt(17,4,15,3,-1,-2)
        expected_G = -2817*np.sqrt(2945)/(731786*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value6(self):
        G = am.gaunt(20,20,40,-1,1,0)
        expected_G = 28384503878959800/(74029560764440771*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value7(self):
        G = am.gaunt(20,15,35,3,2,-5)
        expected_G = -10030817250*np.sqrt(450861502)/(526306917173933*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-16)

    def test_gaunt_confirm_value8(self):
        G = am.gaunt(29,29,34,10,-5,-5)
        expected_G = 1821867940156*np.sqrt(22134)/(215552371055153321*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-12)

    def test_gaunt_confirm_value9(self):
        G = am.gaunt(30,25,35,29,-17,-12)
        expected_G = 39729855*np.sqrt(29557990690363590)/(270587018984128637*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-12)

    def test_gaunt_confirm_value10(self):
        G = am.gaunt(32,13,45,4,-10,6)
        expected_G = 24605*np.sqrt(3127920717610)/(5876912163037*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-12)

    def test_gaunt_confirm_value11(self):
        G = am.gaunt(50,45,25,0,0,0)
        expected_G = 188620746087887071920*np.sqrt(468741)/(1371670599417077008727051*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-11)

    def test_gaunt_confirm_value12(self):
        G = am.gaunt(60,58,118,3,-2,-1)
        expected_G = -6005349369483935328764225756762765828642177466132*np.sqrt(37009780249)/(3204442318367261236376864832217203967591895154546699575*np.sqrt(np.pi))
        assert G == pytest.approx(expected_G, abs=1.0E-4)

    def test_gaunt_confirm_value13(self):
        G = am.gaunt(84,93,53,32,-14,-18)
        try:
            expected_G = 15949506183362618302382939036*np.sqrt(3505135494518220575103)/(17709297921203169674698082092504932688601*np.sqrt(np.pi))
        except TypeError:
            from math import sqrt, pi
            # numpy sqrt method does not support very large integers
            expected_G = 15949506183362618302382939036*sqrt(3505135494518220575103)/(17709297921203169674698082092504932688601*sqrt(pi))
        assert G == pytest.approx(expected_G, abs=1.0E-6)

    # Type tests

    def test_gaunt_returns_float_type(self):
        G = am.gaunt(17,4,15,3,-1,-2)
        self.assertIsInstance(G, float)

    def test_gaunt_raises_TypeError_when_sum_m_is_not_zero(self):
        # Rule: m1+m2+m3=0
        # These couplings: m1+m2+m3 != 0
        with self.assertRaises(TypeError):
            am.gaunt(17,4,15,-3,-1,-2)
        with self.assertRaises(TypeError):
            am.gaunt(17,4,15,1,2,3)

    def test_gaunt_is_approx_zero_when_sum_l_is_not_even(self):
        # Rule: (l1+l2+l3)%2==0 (i.e. (l1+l2+l3) mod 2 = 0)
        # These couplings: (l1+l2+l3)%2==1
        G = am.gaunt(17,4,14,3,-1,-2)
        assert G == pytest.approx(0)
        G = am.gaunt(10,9,12,9,3,-12)
        assert G == pytest.approx(0)

    def test_gaunt_raises_TypeError_when_triangle_rule_is_not_satisfied(self):
        # Triangle condition: |l2-l3| <= l1 <= l2+l3
        with self.assertRaises(TypeError):
            # Violation: |l2-l3| > l1
            am.gaunt(1,10,12,9,3,-12)
        with self.assertRaises(TypeError):
            # Violation: l1 > l2+l3
            am.gaunt(24,10,12,9,3,-12)

    def test_gaunt_raises_TypeError_with_wrong_number_args(self):
        # No args
        with self.assertRaises(TypeError):
            am.gaunt()
        # Too few args
        with self.assertRaises(TypeError):
            am.gaunt(17,4,15,3,-1)
        # Too many args
        with self.assertRaises(TypeError):
            am.gaunt(17,4,15,3,-1,-2,-8)
            
