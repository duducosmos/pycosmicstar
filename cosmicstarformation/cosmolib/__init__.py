#!/usr/bin/env python

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"

"""
DISCLAIMER:

A FORTRAN wrapper library for cosmology analisys in Python

cosmolib is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.
cosmolib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


AUTHOR:

    Eduardo S. Pereira
    email: pereira.somoza@gmail.com


Initialization:

   import cosmolib
   omegab = 0.04
   omegam = 0.24
   omegal = 0.7
   h = 0.7
   myUniverse = cosmolib.init(omegab,omegam,omegal,h)

FUNCTIONS:
    dtdz(z) : Time and Redshift Relation for S{Lambda}CDM Universe
    dtdzCG : Time and Redshift Relation for  Chaplygin Gas
    rz(z): comove distancy  for S{Lambda}CDM Universe
    rzGC(z) : comove distancy  for  Chaplygin Gas
    dr_dz(z) : variation of comove distance with redshift for
    S{Lambda}CDM Universe
    drGC_dz(z) : variation of comove distance with redshift  for  Chaplygin Gas
    dV_dz(z) : variation of comove volume with redshift for S{Lambda}CDM
    Universe
    age(z):  Age of the Universe for S{Lambda}CDM Universe
    ageCG(z) : Age of the Universe for  Chaplygin Gas

    sigma: the variance of the linear density field.
    dsigma2_dk: Derivative of the variance of the linear density field
                         with respect to the scala factor
    grow: Growth function.

    rhodm : Evolution of the dark matter density
    rhobr : Evolution of the barionic matter density.

    NUMERICAL METHODS:
        rk4_in:
            4th-order Runge-Kutta method for solving the
            initial value problem { y}' = { F(x,{ y} )} , where
            { y} = { y[0],y[1],...y[n-1]} .
            ARGUMENTS:
                y: nitial conditions.
                Xarray : Array with the x value for all range
                Yarray: Output Array with the solutions of Y'
                fun: user-supplied function that returns the
                                 array F(x,y) = { y'[0],y'[1],...,y'[n-1]} .
            RETURN:
                The integrated numerical function
        romberg:
            Romberg Integration
            ARGUMENTS:
                func : Function to be integrated
                a : start point
                b : end point
                tol : tolerance
            RETURN:
                The integrated value

        locate:
            Localiza a posicao de dado ponto a partir de dois adjacentes.
            ARGUMENTS:
                func --- function or entry table
                xx   --- entry table
                n    --- n point in the table
                x    ---  value in x that is related to y
            RETURN
                j    ---  x,y position

        dfridr:
            Gives the derivate of func with respect to x
            Arguments:
                func - Function to be derived
                x - point in x where the derivative is analysed
                h - step for diferential
                error - internal parameter for function error
            RETURN:
                The derived value of the funtion in the x point


"""

__all__ = ['lcdmlib', 'chaplyginlib', 'numericallib']
