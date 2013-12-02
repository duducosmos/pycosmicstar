# -*- coding: utf-8 -*-
## module run_kut4

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"


"""4th-order Runge-Kutta method for solving the initial value problem
X,Y = integrate(F,x,y,xStop,h).
   4th-order Runge-Kutta method for solving the
   initial value problem { y} ’ = { F(x,{ y} )} , where
    { y} = { y[0],y[1],...y[n-1]} .
    x,y    = initial conditions.
    xStop = terminal value of x.
    h      = increment of x used in integration.
    F      = user-supplied function that returns the
              array F(x,y) = { y’[0],y’[1],...,y’[n-1]} .


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

from numpy import array


def rk4_int(F, x, y, xStop, h):

    def run_kut4(F, x, y, h):
        # Computes increment of y from Eqs. (7.10)
        K0 = h * F(x, y)
        K1 = h * F(x + h / 2.0, y + K0 / 2.0)
        K2 = h * F(x + h / 2.0, y + K1 / 2.0)
        K3 = h * F(x + h, y + K2)
        #print 'run',K0,K1,K2,K3
        return (K0 + 2.0 * K1 + 2.0 * K2 + K3) / 6.0

    X = []
    Y = []
    X.append(x)
    Y.append(y)

    while x < xStop:
        h = min(h, xStop - x)
        y = y + run_kut4(F, x, y, h)
        x = x + h

        X.append(x)
        Y.append(y)

    return array(X), array(Y)
