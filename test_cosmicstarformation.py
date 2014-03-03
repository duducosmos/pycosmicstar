#!/usr/bin/env python3
# *-* Coding: UTF-8 *-*

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"


"""Cosmic Star Formation Rate

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

from pycosmicstar.cosmicstarformation import cosmicstarformation
from pycosmicstar.lcdmcosmology import lcdmcosmology


class test_cosmicstarformation(unittest.TestCase):

    myCosmicStar = cosmicstarformation(cosmology=lcdmcosmology, tau=2.5)

    def test_cosmicStarsDensity(self):
        self.assertEqual(
            round(self.myCosmicStar.cosmicStarFormationRate(4.5), 3), 0.159)

    def test_gasDensityInStructures(self):
        self.assertEqual(
            round(self.myCosmicStar.gasDensityInStructures(4.5)[0], 2),
            397760539.19)

    def test_phi(self):
        self.assertEqual(round(self.myCosmicStar.phi(1e3), 11), 1.513e-08)

if(__name__ == "__main__"):
    unittest.main()
