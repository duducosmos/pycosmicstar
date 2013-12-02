#!/usr/bin/env python
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

from structures import structures
from lcdmcosmology import lcdmcosmology


class test_structures(unittest.TestCase):

    myStructures = structures(lcdmcosmology)

    def test_funcMassST(self):
        self.assertEquals(round(self.myStructures.funcMassST(9.0, 1.0), 11),
                          8.45e-09)

    def test_fstm(self):
        self.assertEquals(round(self.myStructures.fstm(6.0), 2),
                          515.94)

    def test_halos_n(self):
        self.assertEquals(round(self.myStructures.halos_n(0.0), 2),
                          23581005304.07)

    def test_fbstruc(self):
        self.assertEquals(round(self.myStructures.fbstruc(0.0), 2),
                          0.73)

    def test_numerical_density_halos(self):
        self.assertEquals(
            round(self.myStructures.numerical_density_halos(0.0), 7),
                          6.76e-05)

    def test_abt(self):
        self.assertEquals(round(self.myStructures.abt(1.0), 4),
                            0.0082)

    def test_creatCachDiretory(self):
        self.assertTrue(self.myStructures.getCacheDir()[0],
        "The directory not Exist")


if(__name__ == "__main__"):
    unittest.main()