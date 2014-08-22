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

"""Abstract Class of cosmological models.

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


class Cosmology:

    def dt_dz(self, z):
        """Return the relation between the cosmic time and the redshift
        """
        raise NotImplementedError('I need to be implemented!')

    def dr_dz(self, z):
        """Return the comove-radii"""
        raise NotImplementedError('I need to be implemented!')

    def dV_dz(self, z):
        """Return the comove volume"""
        raise NotImplementedError('I need to be implemented!')

    def rodm(self, z):
        """Return the Dark Matter Density"""
        raise NotImplementedError('I need to be implemented!')

    def robr(self, z):
        """Return the Barionic Matter Density"""
        raise NotImplementedError('I need to be implemented!')

    def H(self, z):
        """Return the Hubble Parameter"""
        raise NotImplementedError('I need to be implemented!')

    def dgrowth_dt(self, z):
        """Return the derivative of growth function of the
         primordial perturbations"""
        raise NotImplementedError('I need to be implemented!')

    def growthFunction(self, z):
        """Return the growth function of the primordial perturbations"""
        raise NotImplementedError('I need to be implemented!')

    def dsigma2_dk(self, kl):
        """"Return the integrating of sigma(M,z) for a top-hat filtering.
        In z = 0 return sigma_8, for z > 0 return sigma(M,z)
        """
        raise NotImplementedError('I need to be implemented!')

    def sigma(self):
        """Return  the variance of the linear density field.
        As pointed out by Jenkis et al. (2001),
        this definition of the mass function has the advantage that it does
        not explicitly depend on redshift, power spectrum or cosmology.
        """
        raise NotImplementedError('I need to be implemented!')

    def age(self, z):
        """Return the age of the Universe for a given z
        """
        raise NotImplementedError('I need to be implemented!')

    def omegamz(self, z):
        raise NotImplementedError('I need to be implemented!')

    def setCosmologicalParameter(self):
        raise NotImplementedError('I need to be implemented!')

    def getCosmologicalParameter(self):
        raise NotImplementedError('I need to be implemented!')

    def getTilt(self):
        raise NotImplementedError('I need to be implemented!')

    def getRobr0(self):
        raise NotImplementedError('I need to be implemented!')

    def getRodm0(self):
        raise NotImplementedError('I need to be implemented!')
