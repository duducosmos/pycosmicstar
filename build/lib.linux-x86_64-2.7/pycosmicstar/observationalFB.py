# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# *-* Coding: UTF-8 *-*

from __future__ import division, absolute_import

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"


"""Observational Cosmic Star Formation Rate

This file is part of pystar.
copyright : Eduardo dos Santos Pereira

pystar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.
pystar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

"""

from numpy import array, loadtxt, log10


class ObservationalFB:

    def __init__(self):
        arq = open("./data/bonamente_2008.dat", 'r')
        self.data = loadtxt(arq, delimiter=",")

    def fb(self):
        """Return the redshift and the fraction of gas  from observational data
        """
        return self.data[:, 0], self.data[:, 4]

    def fbError(self):
        """
        Return the asymetric error in the fraction of gas
        """
        yerr = array([self.data[:, 5], self.data[:, 6]])
        return yerr

    def massTotError(self):
        yerr = array([self.data[:, 2], self.data[:, 3]])
        return yerr

    def minMaxMass(self):
        p = self.data[:, 1] + self.data[:, 2]
        n = self.data[:, 1] + self.data[:, 3]

        ma = p.max()
        mi = n.min()

        return [log10(ma * 1e14), log10(mi * 1e14)]