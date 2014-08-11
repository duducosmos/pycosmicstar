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


"""Cosmic Star Formation Rate

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

from numpy import log10, sqrt, array
from numpy.numarray import zeros, Float64
from .structures import structures
import scipy.interpolate as spint
from scipy.integrate import romberg
from .run_kut4 import rk4_int

import sys
pyversion = sys.version_info
if pyversion[0] >= 3:
    from . import filedict
else:
    print("Importing filedict for python2.7")
    from . import filedict_old as filedict


class cosmicstarformation(structures):
    """The Cosmic Star Formation rate
    The model used to develop this class was presented by the first time
    in the article of Pereira and Miranda (2010) - (MNRAS, 401, 1924, 2010).

    The cosmologic background model is passed as a instance parameter:
        cosmology

    Keyword arguments:
        tau -- (default - 2.5) time scale, in Gyr, of the CSFR.

        eimf -- (default 1.35) exponent of the Initial Mass Function

        nsch -- (default 1) the normalization factor in the CSFR model

        imfType -- (default S - Salpeter) the Initial Mass Function Type.
                   Possible values:
                       S: Salpeter
                       K: Kroupa

        lmin -- (default 6.0) log10 of the minal mass of the dark halo
                            where it is possible to have star formation.

        zmax -- (defaul 20.0) - the maximum redshift to be considered

        omegam -- (default 0.24) - The dark matter parameter

        omegab -- (default 0.04) - The barionic parameter

        omegal -- (default 0.73) - The dark energy parameter

        h -- (default 0.7) - The h of the Hubble constant (H = h * 100)

        massFunctionType -- (default \"ST\") The type of mass
        function of dark matter halos used. Possibles values:
             \"ST\" for Seth and Thormen mass function.
             \"TK\" for Tinker et al. mass function.
    """

    def __init__(self, cosmology,
                       tau=2.29, eimf=1.35, nsch=1, lmin=6.0, zmax=20.0,
                       imfType="S", **kwargs):

        structures.__init__(self, cosmology, **kwargs)

        cacheFile = self._cacheFIle + "_CSFR_" + str(tau)\
                 + "_" + str(eimf) + "_" + str(nsch) + "_" + imfType

        self._cache_dictS = filedict.FileDict(filename=cacheFile + ".cache")

        tau = tau * 1.0e9
        self._cc = self._tck_ab[1]

        #Cosmic Star Formation Rate normalization
        self.__esnor = 1.0
        self.__nsch = nsch
        self.__tau = tau

        self.imfType = imfType

        self.__aminf1 = 2.5e1
        self.__amsup1 = 1.4e2
        self.__amin = 1.0e-1

        self.__imfDict = {"S": self.__imfSalpeter,
                          "K": self.__imfKroupa
                          }

        self.__anorm1 = None
        self.__eimf = eimf
        self.__eimf0 = eimf - 1.0

        try:

            self.__astar = self._cache_dictS['astar']
            self.__csfr = self._cache_dictS['csfr']
            self.__rho_gas = self._cache_dictS['rho_gas']
            print("Data CSFR in Cache")

        except:
            self.__csfr, self.__rho_gas, self.__astar = self.__sfr()

        tck_sf = spint.splrep(self.__astar, self.__csfr)
        self.__cs = tck_sf[1]

        tck_sg = spint.splrep(self.__astar, self.__rho_gas)
        self.__cs2 = tck_sg[1]

    def getIMFDict(self):
        """
        Return a list with keys and functions of IMF's
        """
        dic = []
        for key, value in list(self.__imfDict.items()):
            dic.append([key, value])
        return dic

    def putIMFDict(self, key, value):
        """
        Put a new term in the imf Dictionary
        """

        self.__imfDict[key] = value

    def __spn(self, a):
        """Return the interpolated, by cubi-spline, barionic accretion
        rate into structures
        """
        j = -1

        while(1):
            j = j + 1
            if(a < self._ascale[j + 1] and
               a >= self._ascale[j] and
               j <= (len(self._ascale) - 3)
               ):
                resp = (self._cc[j] / 6.0) * \
                       (((a - self._ascale[j + 1]) ** 3.0) /
                       (self._ascale[j] - self._ascale[j + 1])
                       - (a - self._ascale[j + 1]) *
                       (self._ascale[j] - self._ascale[j + 1])) \
                       - (self._cc[j + 1] / 6.0) * \
                       (((a - self._ascale[j]) ** 3.0) /
                       (self._ascale[j] - self._ascale[j + 1])
                       - (a - self._ascale[j]) *
                       (self._ascale[j] - self._ascale[j + 1])) \
                       + (self._abt2[j] * (a - self._ascale[j + 1])
                       - self._abt2[j + 1] * (a - self._ascale[j]))\
                       / (self._ascale[j] - self._ascale[j + 1])
                return resp

            elif(j >= (len(self._ascale) - 2) and j <= len(self._ascale)):
                return self._abt2[j]
            elif(a < self._ascale[0]):
                raise NameError("Error in the spline function")
                break

    def __fcn(self, a, rho_g):
        """Return the numerical function to be integrated to calculate
        the mass density of barions into structures.
        """

        z = 1.0 / a - 1.0

        if(z < 0.0):
            z = 0.0

        tage = self._cosmology.age(z)

        age01 = 4.0 * log10(tage) - 2.704e+01
        age02 = (3.6 - sqrt(age01)) / 2.0

        mi_1 = 1.0e+01 ** age02
        yr = self.__massEjected(mi_1)

        if(self.__nsch == 1.0):
            sexp = (1.0 - yr) / self.__tau
        else:
            sexp = (1.0 - yr) / self.__tau /\
             self._cosmology.getRobr0() ** (self.__nsch - 1.0)

        F = zeros(1, type=Float64)
        F[0] = (-sexp * (rho_g[0]) ** self.__nsch
                + self.__esnor * self.__spn(a) /
                self._cosmology.getRobr0()) \
                * self._cosmology.dt_dz(z) / a ** 2.0
        return F

    def __csfr_gas(self, rg):
        """Return the Cosmic Star Formation Rate
        from the barionic gas into structures
        """
        if(self.__nsch == 1):
            roes = rg / self.__tau
        else:
            roes = (rg ** self.__nsch) / self.__tau\
            / self._cosmology.getRobr0() ** (self.__nsch - 1.0)
        return roes

    def __sfr(self):
        """Return the Cosmic Star Formation rate, density of barionic
        gas into structures
        """

        #Normalization of the cosmic star formation rate
        rho_g0 = array([1.0e-9])
        a0 = self._ascale[0]
        nf = len(self._ascale) - 1
        af = self._ascale[nf]
        step = (af - a0) / 100.0

        A, R_g = rk4_int(self.__fcn, a0, rho_g0, af, step)

        ng = len(A) - 1

        if(self.__nsch == 1):
            roes = R_g[ng] / self.__tau
        else:
            roes = (R_g[ng] ** self.__nsch) / self.__tau /\
            self._cosmology.getRobr0() ** (self.__nsch - 1.0)

        self.__esnor = 1.62593696e-2 / roes

        rho_g0 = array([1.0e-9])
        a0 = 1.0 / (self._zmax + 1.0)
        nf = len(self._ascale) - 1
        af = self._ascale[nf]
        step = (af - a0) / 5000.
        A, R_g = rk4_int(self.__fcn, a0, rho_g0, af, step)

        rho_s = self.__csfr_gas(R_g)

        self._cache_dictS['astar'] = A
        self._cache_dictS['csfr'] = rho_s
        self._cache_dictS['rho_gas'] = R_g

        return rho_s, R_g, A

    def __massEjectedSalpeter(self, m_min):
        if(self.__anorm1 is None):
            self.__imfSalpeter(10)
        amexp1 = (1.0 / m_min) ** self.__eimf0
        amexp2 = (1.0 / self.__amsup1) ** self.__eimf0
        amexp3 = (1.0 / 8.0) ** self.__eimf0
        amexp4 = (1.0 / self.__aminf1) ** self.__eimf0
        amexp5 = (1.0 / m_min) ** self.__eimf
        amexp6 = (1.0 / 8.0) ** self.__eimf
        amexp7 = (1.0 / 10.0) ** self.__eimf
        amexp8 = (1.0 / self.__aminf1) ** self.__eimf
        amexp9 = (1.0 / self.__amsup1) ** self.__eimf

        yrem1 = (amexp1 - amexp2) / self.__eimf0
        yrem2 = 1.156e-01 * (amexp1 - amexp3) / self.__eimf0
        yrem3 = 1.3e+01 * (amexp4 - amexp2) / self.__eimf0 / 2.4e+01
        yrem4 = 4.551e-01 * (amexp5 - amexp6) / self.__eimf
        yrem5 = 1.35e+00 * (amexp6 - amexp7) / self.__eimf
        yrem6 = 1.40e+00 * (amexp7 - amexp8) / self.__eimf
        yrem7 = 6.5e+01 * (amexp8 - amexp9) / self.__eimf / 6.0

        yr = self.__anorm1 * (yrem1 - yrem2 - yrem3
                            - yrem4 - yrem5 - yrem6 + yrem7)

        return yr

    def __massEjected(self, m_min):
        """
        Return the mass integration of the mass ejected by the collapse of the
        star
        """
        if(self.imfType == "S"):
            return self.__massEjectedSalpeter(m_min)

        else:
            mEject = (romberg(self.__mPhi, m_min, self.__amsup1, tol=1.48e-04)
                 - romberg(self.__mrPhi, m_min, self.__amsup1, tol=1.48e-04)
                      )
        return mEject

    def __mPhi(self, m):
        return m * self.phi(m)

    def __mrPhi(self, m):
        return self.remnant(m) * self.phi(m)

    def remnant(self, m):
        """
        Return the remnant mass of the object after the colapse of the star
        with mass m
        """

        if(m > 0 and m < 1):
            return 0
        elif(m >= 1 and m <= 8):
            return 0.1156 * m + 0.4551
        elif(m > 8 and m <= 10):
            return 1.35
        elif(m > 10 and m < 25):
            return 1.4
        elif(m >= 25 and m <= 145):
            return (13.0 / 24.0) * (m - 20)
        else:
            raise NameError("Error: Out of the mass range...")

    def __imfKroupa(self, m):
        if(self.__anorm1 is None):
            alpha0 = 0.3
            alpha1 = 1.3
            alpha2 = 2.3
            alpha3 = alpha2

            k0 = 1
            k1 = k0 * 0.08
            k2 = k1 * 0.5
            k3 = k2
            A = [k0 * (self.__amsup1 ** (1.0 - alpha0)
                     - self.__amin ** (1.0 - alpha0)) / (1.0 - alpha0),
                k1 * (self.__amsup1 ** (1.0 - alpha1)
                     - self.__amin ** (1.0 - alpha1)) / (1.0 - alpha1),
                k2 * (self.__amsup1 ** (1.0 - alpha2)
                     - self.__amin ** (1.0 - alpha2)) / (1.0 - alpha2),
                k3 * (self.__amsup1 ** (1.0 - alpha3)
                     - self.__amin ** (1.0 - alpha3)) / (1.0 - alpha3)
                ]

            self.__anorm1 = 1 / sum(A)

        if(m > 0.08 and m <= 0.5):
            return self.__anorm1 * m ** (-1.3)
        elif(m > 0.5):
            return self.__anorm1 * m ** (-2.3)
        else:
            raise NameError("Mass out of the range")

    def __imfSalpeter(self, m):

        if(self.__anorm1 is None):
            self.__amin = 1.0e-1
            self.__anorm1 = self.__eimf0 / (1.0 / self.__amin ** self.__eimf0
                             - 1.0 / self.__amsup1 ** self.__eimf0)

        return self.__anorm1 * m ** (-(1.0 + self.__eimf))

    def phi(self, m):
        """Return the Initial Mass Function
        """
        try:
            return self.__imfDict[self.imfType](m)
        except:
            raise NameError("No Defined Initial Mass Function")

    def cosmicStarFormationRate(self, z):
        """Return the Cosmic Star Formation rate as a function of z
        """

        a = 1.0 / (1.0 + z)

        j = -1

        while(1):
            j = j + 1
            if(a < self.__astar[j + 1] and
               a >= self.__astar[j] and
               j <= (len(self.__astar) - 3)):
                resp = (self.__cs[j] / 6.0) * \
                       (((a - self.__astar[j + 1]) ** 3.0) /
                       (self.__astar[j] - self.__astar[j + 1])
                       - (a - self.__astar[j + 1]) *
                       (self.__astar[j] - self.__astar[j + 1])) \
                       - (self.__cs[j + 1] / 6.0) * \
                       (((a - self.__astar[j]) ** 3.0) /
                       (self.__astar[j] - self.__astar[j + 1])
                       - (a - self.__astar[j]) *
                       (self.__astar[j] - self.__astar[j + 1])) \
                       + (self.__csfr[j] * (a - self.__astar[j + 1])
                       - self.__csfr[j + 1] * (a - self.__astar[j])) / \
                        (self.__astar[j] - self.__astar[j + 1])
                return resp[0]

            elif(j >= (len(self.__astar) - 2) and j <= len(self.__astar)):
                return self.__csfr[j]
            elif(a < self.__astar[0]):
                raise NameError("Error in spline csfr")
                break

    def gasDensityInStructures(self, z):
        """Return the barionic gas density into structures
        """
        a = 1.0 / (1.0 + z)

        j = -1

        while(1):
            j = j + 1
            if(a < self.__astar[j + 1] and
               a >= self.__astar[j] and
               j <= (len(self.__astar) - 3)):
                resp = (self.__cs2[j] / 6.0) * ((
                       (a - self.__astar[j + 1]) ** 3.0) /
                       (self.__astar[j] - self.__astar[j + 1])
                       - (a - self.__astar[j + 1]) *
                       (self.__astar[j] - self.__astar[j + 1])) \
                       - (self.__cs2[j + 1] / 6.0) * \
                       (((a - self.__astar[j]) ** 3.0) /
                       (self.__astar[j] - self.__astar[j + 1])
                       - (a - self.__astar[j]) *
                       (self.__astar[j] - self.__astar[j + 1])) \
                       + (self.__rho_gas[j] * (a - self.__astar[j + 1])
                       - self.__rho_gas[j + 1] * (a - self.__astar[j])) \
                       / (self.__astar[j] - self.__astar[j + 1])
                return resp

            elif(j >= (len(self.__astar) - 2) and j <= len(self.__astar)):
                return self.__rho_gas[j]
            elif(a < self.__astar[0]):
                raise NameError("Error spline gas density")
                break
