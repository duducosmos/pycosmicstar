#!/usr/bin/env python3
# *-* Coding: UTF-8 *-*
__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"

"""
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

from pycosmicstar.structures import Structures
from pycosmicstar.lcdmcosmology import Lcdmcosmology


class test_structures(unittest.TestCase):

    myStructures = Structures(cosmology=Lcdmcosmology)

    def test_massFunction(self):
        self.assertEqual(round(self.myStructures.massFunction(9.0, 1.0), 11),
                          8.45e-09)

    def test_fstm(self):
        self.assertEqual(round(self.myStructures.fstm(6.0), 2),
                          515.94)

    def test_halos_n(self):
        self.assertEqual(round(self.myStructures.halos_n(0.0), 2),
                          23584651682.25)

    def test_fbstruc(self):
        self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          0.73)

    def test_numerical_density_halos(self):
        self.assertEqual(
            round(self.myStructures.numerical_density_halos(0.0), 7),
                          6.76e-05)

    def test_abt(self):
        self.assertEqual(round(self.myStructures.abt(1.0), 4),
                            0.0082)

    def test_creatCachDiretory(self):
        self.assertTrue(self.myStructures.getCacheDir()[0],
        "The directory not Exist")

    def test_setDeltaHTinker(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="TK")
        self.assertTrue(self.myStructures.setDeltaHTinker(200))

    def test_massfunctioR(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="R")
        self.assertEqual(round(self.myStructures.massFunction(9.0, 1.0), 11),
                          8.45e-09)

    def test_fbstrucR(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="R")
        self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          0.71)

    def test_massfunctioPS(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="PS")
        self.assertEqual(round(self.myStructures.massFunction(9.0, 1.0), 11),
                          4.54e-09)

    def test_fbstrucPS(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="PS")
        self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          0.82)

    def test_massfunctioWT(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="WT1")
        self.assertEqual(round(self.myStructures.massFunction(9.0, 1.0), 9),
                          2.60e-08)

    def test_fbstrucWT(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="WT1")
        self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          1.15)

    def test_massfunctioWT2(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="WT2")
        self.assertEqual(round(self.myStructures.massFunction(9.0, 1.0), 10),
                          1.76e-08)

    def test_fbstrucWT2(self):
        self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       massFunctionType="WT2")
        self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          0.87)

    #def test_massfunctioJK(self):
        #self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       #massFunctionType="JK")
        #self.assertEqual(self.myStructures.massFunction(9.0, 1.0),
                          #8.45e-09)

    #def test_fbstrucJK(self):
        #self.myStructures = Structures(cosmology=Lcdmcosmology,
                                       #massFunctionType="JK")
        #self.assertEqual(round(self.myStructures.fbstruc(0.0), 2),
                          #0.47)

if(__name__ == "__main__"):
    unittest.main()
