#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"

"""Unit Test module for the cosmological model

This file is part of cosmicstar.
copyright : Eduardo dos Santos Pereira

cosmicstar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.
cosmicstar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

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
