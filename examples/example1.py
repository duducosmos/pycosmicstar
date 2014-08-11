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


"""Example 1.

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

import matplotlib.pyplot as plt
from numpy import arange, log
from pycosmicstar.lcdmcosmology import lcdmcosmology

lcdmUniverser = lcdmcosmology(omegam=0.24,
                              omegab=0.04,
                              omegal=0.73,
                               h=0.7)


z = arange(0, 10.5, 0.1)

plt.plot(z, [lcdmUniverser.dt_dz(zi) for zi in z])
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"$\frac{dt}{dz}$(yr)")
plt.show()


plt.plot(z, [lcdmUniverser.age(zi) for zi in z])
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"$t$ (yr)")
plt.show()


plt.plot(z, [lcdmUniverser.H(zi) for zi in z])
plt.xlabel(r"$z$ - Redshift")
plt.ylabel(r"$H(z)$ ")
plt.show()
