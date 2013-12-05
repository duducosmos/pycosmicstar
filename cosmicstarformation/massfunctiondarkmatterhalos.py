#!/usr/bin/env python3.3
# *-* Coding: UTF-8 *-*

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"


"""Mass Function of dark matter halos.
A colection of mass function of dark matter halos

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

from numpy import sqrt, exp, pi


def pressSchechterMassFunction(nu):
    """Return the value of Press-Schechter (1974) mass function.
    Keyword arguments:
        nu - S{delta}_{c}/S{sigma}
    """

    fmass = sqrt(2.0 / pi) * (nu) * exp(-0.5 * (nu) ** 2)
    return fmass