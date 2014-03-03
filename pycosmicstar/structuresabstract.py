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

"""Abstract Class of like Press-Schechter formalism

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


class structuresabstract:

    def massFunction(self, lm, z):
        """Return the mass function of dark halos
        """
        raise NotImplementedError('I need to be implemented!')

    def fstm(self, lm):
        '''Numerical function that return the value of sigm that
        will be used by dfridr to calculate d_sigma_dlog10(m).'''
        raise NotImplementedError('I need to be implemented!')

    def halos_n(self, z):
        """Return the integral of the mass function of dark halos multiplied
        by mass in the range of log(M_min) a log(M_max)
        """
        raise NotImplementedError('I need to be implemented!')

    def fbstruc(self, z):
        """Return the faction of barions into structures
        """
        raise NotImplementedError('I need to be implemented!')

    def numerical_density_halos(self, z):
        """Return the numerial density of dark halos
        within the comove volume
        """
        raise NotImplementedError('I need to be implemented!')

    def abt(self, a):
        """Return the accretion rate of barionic matter into strutures
        """
        raise NotImplementedError('I need to be implemented!')

    def cosmicStarFormationRate(self, z):
        raise NotImplementedError("I need to be implemented!")

    def gasDensityInStructures(self, z):
        raise NotImplementedError("I need to be implemented!")

    def getCacheDir(self):
        """Return True if the cache directory existe and false else.
        """
        raise NotImplementedError('I need to be implemented!')
