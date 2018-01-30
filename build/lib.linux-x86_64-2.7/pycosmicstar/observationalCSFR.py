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

from numpy import array, loadtxt


class ObservationalCSFR:

    def __init__(self):
        arq = open("./data/hopkins_2004.dat", 'r')
        self.data = loadtxt(arq, delimiter=",")

    def csfredshift(self):
        """Return the redshift and the CSFR from
        observational data
        """
        return self.data[:, 0], self.data[:, 1]

    def errorData(self):
        """Return the asymetric errors in the redshif and CSFR
        respectively
        """
        xerr = array([self.data[:, 0] - self.data[:, 2],
                self.data[:, 3] - self.data[:, 0]])
        yerr = array([self.data[:, 1] - self.data[:, 4],
                      self.data[:, 5] - self.data[:, 1]])
        return xerr, yerr
