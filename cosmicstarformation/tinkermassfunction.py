#!/usr/bin/env python
# *-* Coding: UTF-8 *-*
from __future__ import print_function

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"

""" Tinker mass function (Tinker et al. 2008)

This function was adapted from the work of:
    S.G. Murray et al. 2013. Astronomy and Computing. 3-4. 23-34.

source of the original (https://github.com/steven-murray/hmf)

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


from numpy import array, exp, log, logical_or, nan
from scipy.interpolate import InterpolatedUnivariateSpline as spline


def funcMass(self):
    """Return the mass function of dark halos of
Tinker mass function (Tinker et al. 2008)

This function was adapted from the work of:
    S.G. Murray et al. 2013. Astronomy and Computing. 3-4. 23-34.
    source of the original (https://github.com/steven-murray/hmf)

Keyword arguments:
    lm -- log10 of the mass of the dark halo
    z -- redshift
    """

    #The Tinker function is a bit tricky - we use the code from
    #http://cosmo.nyu.edu/~tinker/massfunction/MF_code.tar
    #to aid us.
    delta_virs = array([200, 300, 400, 600, 800, 1200, 1600, 2400, 3200])
    A_array = array([1.858659e-01,
                        1.995973e-01,
                        2.115659e-01,
                        2.184113e-01,
                        2.480968e-01,
                        2.546053e-01,
                        2.600000e-01,
                        2.600000e-01,
                        2.600000e-01])

    a_array = array([1.466904e+00,
                        1.521782e+00,
                        1.559186e+00,
                        1.614585e+00,
                        1.869936e+00,
                        2.128056e+00,
                        2.301275e+00,
                        2.529241e+00,
                        2.661983e+00])

    b_array = array([2.571104e+00,
                        2.254217e+00,
                        2.048674e+00,
                        1.869559e+00,
                        1.588649e+00,
                        1.507134e+00,
                        1.464374e+00,
                        1.436827e+00,
                        1.405210e+00])

    c_array = array([1.193958e+00,
                        1.270316e+00,
                        1.335191e+00,
                        1.446266e+00,
                        1.581345e+00,
                        1.795050e+00,
                        1.965613e+00,
                        2.237466e+00,
                        2.439729e+00])
    A_func = spline(delta_virs, A_array)
    a_func = spline(delta_virs, a_array)
    b_func = spline(delta_virs, b_array)
    c_func = spline(delta_virs, c_array)

    A_0 = A_func(self.delta_halo)
    a_0 = a_func(self.delta_halo)
    b_0 = b_func(self.delta_halo)
    c_0 = c_func(self.delta_halo)

    A = A_0 * (1 + self.z) ** (-0.14)
    a = a_0 * (1 + self.z) ** (-0.06)
    alpha = exp(-(0.75 / log(self.delta_halo / 75)) ** 1.2)
    b = b_0 * (1 + self.z) ** (-alpha)
    c = c_0

    vfv = A * ((self.sigma / b) ** (-a) + 1) * exp(-c / self.sigma ** 2)

    if self.cut_fit:
        if self.z == 0.0:
            vfv[logical_or(self.lnsigma / log(10)
                < -0.6, self.lnsigma / log(10) > 0.4)] = nan
        else:
            vfv[logical_or(self.lnsigma / log(10)
            < -0.2, self.lnsigma / log(10) > 0.4)] = nan
    return vfv