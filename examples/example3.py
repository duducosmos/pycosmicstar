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


"""Example 3.

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
from pycosmicstar.cosmicstarformation import cosmicstarformation

from pycosmicstar.lcdmcosmology import lcdmcosmology
from pycosmicstar.observationalCSFR import ObservationalCSFR

import matplotlib.pyplot as plt
from numpy import arange, array


z = arange(0, 20, 0.1)

#Cosmic Star Formation Rate using Tinker et al. dark haloes mass function
myCSFR_TK = cosmicstarformation(cosmology=lcdmcosmology,
                                massFunctionType="TK",
                                delta_halo=200)

#Cosmic Star Formation Rate using Press and Schechter dark haloes mass function
myCSFR_PS = cosmicstarformation(cosmology=lcdmcosmology,
                                massFunctionType="PS")

#Cosmic Star Formation Rate using Seth et al. dark haloes mass function
myCSFR_ST = cosmicstarformation(cosmology=lcdmcosmology)


csfrTK = array([myCSFR_TK.cosmicStarFormationRate(zi) for zi in z])

csfrPS = array([myCSFR_PS.cosmicStarFormationRate(zi) for zi in z])
csfrST = array([myCSFR_ST.cosmicStarFormationRate(zi) for zi in z])

obsCSFR = ObservationalCSFR()
x, y = obsCSFR.csfredshift()
xerr, yerr = obsCSFR.errorData()
plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='.')
plt.plot(z, csfrTK, label="TK")
plt.plot(z, csfrST, label="ST")
plt.plot(z, csfrPS, label="PS")
plt.legend(loc=4)
plt.yscale('log')
plt.ylabel(r'$\dot{\rho}_{*}$( M$_{\odot}$Mpc$^{-3}$yr$^{-1}$)')
plt.xlabel(r'$z$')
plt.show()
