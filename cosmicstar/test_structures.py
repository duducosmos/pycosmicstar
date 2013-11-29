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
        self.assertEquals(self.myStructures.funcMassST(9.0, 1.0),
                          8.4508147954749659e-09)

    def test_fstm(self):
        self.assertEquals(self.myStructures.fstm(6.0), 515.85483654505242)

    def test_halos_n(self):
        self.assertEquals(self.myStructures.halos_n(0.0), 23581044522.84875)

    def test_fbstruc(self):
        self.assertEquals(self.myStructures.fbstruc(0.0), 0.7265184278211807)

    def test_numerical_density_halos(self):
        self.assertEquals(self.myStructures.numerical_density_halos(0.0),
                            6.757090440025199e-05)

    def test_abt(self):
        self.assertEquals(self.myStructures.abt(1.0), 0.0078707904676166528)

    def test_creatCachDiretory(self):
        self.assertTrue(self.myStructures.getCacheDir(),
        "The directory not Exist")



if(__name__ == "__main__"):
    unittest.main()