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


"""Example 2.

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
from pycosmicstar.structures import structures

from pycosmicstar.lcdmcosmology import lcdmcosmology

import matplotlib.pyplot as plt
from numpy import arange, log


#Creating a object considering the
#Sheth & Tormen (1999) mass function.
stUniverser = structures(cosmology=lcdmcosmology,
                         omegam=0.24,
                         omegab=0.04,
                         omegal=0.73,
                         h=0.7,
                         massFunctionType="ST")

#Creating a object considering the
#Tinker et al. (2010) mass function.
tkUniverser = structures(cosmology=lcdmcosmology,
                         omegam=0.24,
                         omegab=0.04,
                         omegal=0.73,
                         h=0.7,
                         massFunctionType="TK",
                         delta_halo=200)


z = arange(0, 10.5, 0.1)

plt.plot(z, [stUniverser.fbstruc(zi) for zi in z], label="ST")
plt.plot(z, [tkUniverser.fbstruc(zi) for zi in z], label="TK")
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"$f_{b}$")
plt.legend()
plt.savefig("./examples/fb.png")
plt.show()

plt.plot(z, [stUniverser.abt(1.0 / (1.0 + zi)) for zi in z], label="ST")
plt.plot(z, [tkUniverser.abt(1.0 / (1.0 + zi)) for zi in z], label="TK")
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"$a_{b}$")
plt.legend()
plt.savefig("./examples/ab.png")
plt.show()


plt.plot(z, [stUniverser.numerical_density_halos(zi) for zi in z], label="ST")
plt.plot(z, [tkUniverser.numerical_density_halos(zi) for zi in z], label="TK")
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"numerical density of halos")
plt.legend()
plt.savefig("./examples/numerical_density_halos.png")
plt.show()
