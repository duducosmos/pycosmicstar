#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

"""
Unit Test module for the cosmological model
"""

import unittest
from lcdmcosmology import lcdmcosmology
from numpy import array


class test_lcdmcosmology(unittest.TestCase):

    myUniverse = lcdmcosmology()

    def test_dtdz(self):
        self.assertEqual(self.myUniverse.dt_dz(0.0), 14185836421.438278)

    def test_drdz(self):
        self.assertEqual(self.myUniverse.dr_dz(0.0), 1.3389180199564208e+23)

    def test_H(self):
        self.assertEqual(self.myUniverse.H(0.0), 68.942004612572731)

    def test_dV_dz(self):
        self.assertEqual(self.myUniverse.dV_dz(0.0), 0.0)

    def test_growthFunction(self):
        self.assertEqual(self.myUniverse.growthFunction(0.0),
                        0.16951226599284788)

    def test_dgrowth_dt(self):
        self.assertEqual(self.myUniverse.dgrowth_dt(0.0),
                        5.949048385364545)

    def test_sigma(self):
        self.assertEqual(self.myUniverse.sigma(9.0)[0][0], 0.95424250943932487)

    def test_rodm(self):
        self.assertEqual(self.myUniverse.rodm(0)[0], 32457599999.999996)

    def test_robr(self):
        self.assertEqual(self.myUniverse.robr(0), 5409600000.0)

    def test_setCosmologicalParameter(self):
        self.assertTrue(self.myUniverse.setCosmologicalParameter(omegam=0.24,
                             omegab=0.04, omegal=0.73, h=0.7),
                        "The model need cosmological parameter")

    def test_getCosmologicalParameter(self):
        self.assertEqual(self.myUniverse.getCosmologicalParameter(),
                        (0.04, 0.24, 0.73, 0.7))

    def test_getDeltaC(self):
        self.assertEqual(self.myUniverse.getDeltaC(), 1.686)

    def test_age(self):
        self.assertEqual(self.myUniverse.age(0), 14428495783.321497)

    def test_getTilt(self):
        self.assertEqual(self.myUniverse.getTilt(), 1.92)

    def test_getRobr0(self):
        self.assertEqual(self.myUniverse.getRobr0(), 5409600000.0)
if(__name__ == "__main__"):
    unittest.main()
