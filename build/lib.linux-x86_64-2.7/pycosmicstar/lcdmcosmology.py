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

"""The Cold Dark Matter plus Cosmological constant Module (LCDM)

Na atual versao usamos a normalizacao do WMAP (sem ondas gravitacionais)
a expressao foi adaptada de Eisenstein e Hu (ApJ 511, 5, 1999)
de forma a fornecer sigma_8 = 0,84.
A fracao de massa dos halos e obtida de Sheth e Tormen (MNRAS 308, 119, 1999)
Todos os modelos consideram Omega_Total = Omega_M + Omega_L = 1,0

"Best Fit" do WMAP-3: omega_m = 0,238, omega_b = 0,042, omega_l = 0,762,
h = 0,734, sigma_8 = 0,744
Veja que sigma_8 pelo WMAP e' obtido atraves da recombinacao. Outras
estimativas (p.e. aglomerados de galaxias) fornecem sigma_8 = 0,84.
Conjunto de dados: WMAP-3: omega_m = 0,238, omega_b = 0,042, omega_l = 0,762
h = 0,734, sigma_8 = 0,84
WMAP-1: omega_m = 0,29, omega_b = 0,44, omega_l = 0,71
h = 0,72, sigma_8 = 0,9

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

from .cosmology import Cosmology
from numpy import sqrt, pi, log, log10, exp, sin, cos
from numpy import zeros
from numpy import float64 as Float64

from scipy.integrate import romberg

cosmolibImportStatus = None
try:
    from .cosmolib import lcdmlib
    cosmolibImportStatus = True
    print('lcdmlib imported')
except:
    cosmolibImportStatus = False
    print('lcdmlib not imported, using pure python version of sigma')


class Lcdmcosmology(Cosmology):
    """The Cold Dark Matter (CDM) plus Cosmolocical Constan (Lambda) -  lcdm

    Keyword arguments:
        omegam -- (default 0.24) - The dark matter parameter

        omegab -- (default 0.04) - The barionic parameter

        omegal -- (default 0.73) - The dark energy parameter

        h -- (default 0.7) - The h of the Hubble constant (H = h * 100)
    """

    def __init__(self, omegam=0.24, omegab=0.04, omegal=0.73, h=0.7):
        self.__omegab = omegab
        self.__omegam = omegam
        self.__omegal = omegal
        self.__h = h
        self.__ct3 = 9.78e+09 / h

        #m / s
        self.__speedOfLight = 3.0e+8

        self.__cosmolibImportStatus = cosmolibImportStatus

        if(self.__cosmolibImportStatus is True):
            self.__lcdmlib = lcdmlib
            self.__lcdmlib.init(omegab, omegam, omegal, h)

        if(self.__omegal >= 0.73):
            tilt = 1.92
        elif(omegal >= 0.7 and omegal <= 0.73):
                tilt = 1.915
        else:
            tilt = 1.8

        self.__tilt = tilt

        h2 = h * h
        h2om = h2 * omegam
        h2br = h2 * omegab
        #m3 kg-1 s-2
        self.G = 6.67e-11
        self.__roc0 = 2.76e+11 * h2
        self.__rodm0 = 2.76e+11 * h2om
        self.__robr0 = 2.76e+11 * h2br
        self.__deltac = 1.686
        self.__hsl = h * sqrt(omegal)
        self.__omegalm = omegal / omegam
        self.__s2pi = pi * sqrt(2.0)
        self.__ut = 1.0 / 3.0
        self.__nr = 14000
        self.__ct0 = 4.0 * pi
        self.__ct1 = self.__ct0 * 2.76e+11 / 3.0
        self.__ct2 = self.__ct1 * h2om
        self.__ct3 = 9.78e+09 / h

        #Power Spectrun Normalization
        self.__anorm = 1.94e-5 * \
                    (self.__omegam ** (-0.785 - 0.05 * log(self.__omegam))) \
                    * exp(- 0.95 * (tilt - 1.) - 0.169 * (tilt - 1.) ** 2.0)\
                    / 2.0 / (pi * pi)

        self.__anorm = ((self.__anorm) ** 2.0) * \
                       ((2997.9 / self.__h) ** (3.0 + tilt))

        self.__gama1 = omegab * (1.0 + sqrt(2.0 * h) / omegam)
        self.__gamam = omegam * (h ** 2.0) / (self.__gama1)
        self.__alfa = 6.4 / self.__gamam / h
        self.__beta = 3.0 / self.__gamam / h
        self.__gama = 1.7 / self.__gamam / h

    def dt_dz(self, z):
        dtdz = self.__ct3 / ((1.0 + z) * sqrt(self.__omegal +
                                    self.__omegam * (1.0 + z) ** 3.0))
        return dtdz

    def dr_dz(self, z):

        #Speed of Light km / s
        vl = 3.0e+5

        #H_{0}/h = 1/s
        hub = 3.25e-18

        drdz = (vl / hub / self.__h) / sqrt(self.__omegam * (1.0 + z) ** 3.0
                                         + self.__omegal)
        return drdz

    def H(self, z):
        """Return the Hubble parameter as a function of z.

        Keyword arguments:
            z -- redshift
        """
        return 100.0 * self.__h * sqrt(self.__omegam * (1.0 + z) ** 3.0
                                     + self.__omegal)

    def dV_dz(self, z):
        """Return the comove volume variation.

        Keyword arguments:
            z -- redshift
        """
        rz = romberg(self.dr_dz, 0.0, z, tol=1.48e-09)
        drdz = self.dr_dz(z)
        dVdz = 4.0 * pi * drdz * rz ** 2.0
        return dVdz

    def dgrowth_dt(self, z):
        """Return the derivative of the growth function with
        respect to  time.

        Keyword arguments:
            z -- redshift
        """
        z1 = 1.0 + z
        ascale = 1.0 / z1
        ascale2 = ascale ** 2.0
        ascale3 = ascale ** 3.0
        ascale4 = ascale * ascale3
        ea = self.__omegam * ascale + self.__omegal * ascale4
        omegamz = self.__omegam * ascale / ea
        omegalz = self.__omegal * ascale4 / ea
        dz1 = 1.0 - omegalz + omegamz ** (4.0 / 7.0) + omegamz / 2.0

        Q = 2.5 * omegamz * ascale
        dea_da = self.__omegam + 4.0 * ascale3 * self.__omegal ** 6.0
        domegamz_da = (self.__omegam / ea ** 2.0) * (ea - ea * dea_da)
        domegalz_da = self.__omegal * (4.0 * ascale3 * ea -
                                       ascale4 * dea_da) / (ea ** 2.0)
        dQ_da = 5.0 * (omegamz + ascale * domegamz_da)
        dP_da = 2.0 * (- domegalz_da + (4.0 / 7.0) *
                      domegamz_da / (omegamz ** (3.0 / 7.0))
                      + domegamz_da / 2.0)
        dadz = ascale2
        dgrowthdt = (dadz) * (dz1 * dQ_da - Q * dP_da) / (dz1 ** 2.0)
        return dgrowthdt

    def growthFunction(self, z):
        """Return the growth function

        Keyword arguments:
            z -- redshift
        """
        z1 = 1.0 + z
        ascale = 1.0 / z1
        ascale3 = ascale ** 3.0
        ascale4 = ascale * ascale3
        ea = self.__omegam * ascale + self.__omegal * ascale4
        omegamz = self.__omegam * ascale / ea
        omegalz = self.__omegal * ascale4 / ea
        dz1 = 1.0 - omegalz + omegamz ** (4.0 / 7.0) + omegamz / 2.0
        growth = (2.5 * omegamz * ascale / dz1) / (pi * sqrt(2.0))

        return growth

    def sigma(self, kmass):
        """Return the sigma.

        Keyword arguments:
            kmass -- mass scale
        """

        if(self.__cosmolibImportStatus is not True):
            return self.__sigma(kmass)
        else:
            return self.__lcdmlib.sigma(self.__anorm,
                                    self.__alfa,
                                    self.__beta,
                                    self.__gama,
                                    self.__ct2,
                                    kmass)

    def dsigma2_dk(self, kl):
        """"Return the integrating of sigma(M,z) for a top-hat filtering.
        In z = 0 return sigma_8, for z > 0 return sigma(M,z)
        """
        k = exp(kl)
        x = self.__escala * k
        pk1 = 1.0 + (self.__alfa * k + (self.__beta * k) ** 1.5
                     + (self.__gama * k) ** 2.0) ** 1.13
        pk2 = 1.0 / pk1
        pdmk = pk2 * (k ** 3.0)
        dsigdk = pdmk * (3.0 * (sin(x) - x * cos(x)) / (x ** 3.0)) ** 2.0
        return dsigdk

    def __sigma(self, kmass):

        n = kmass.size
        km = zeros(n, type=Float64)
        sg = zeros(n, type=Float64)

        for i in range(0, n):
            self.__escala = (kmass[i] / self.__ct2) ** (1.0 / 3.0)
            km[i] = log10(kmass[i])

            t0 = log10(1.0e-7 / self.__escala)
            t1 = log10(1.0e-3 / self.__escala)
            t2 = log10(1.0e+0 / self.__escala)
            t3 = log10(10.0e+0 / self.__escala)
            t4 = log10(100.0e+0 / self.__escala)

            sig2_1 = romberg(self.dsigma2_dk, t0, t1, tol=1.48e-09)
            sig2_2 = romberg(self.dsigma2_dk, t1, t2, tol=1.48e-09)
            sig2_3 = romberg(self.dsigma2_dk, t2, t3, tol=1.48e-09)
            sig2_4 = romberg(self.dsigma2_dk, t3, t4, tol=1.48e-09)

            sg[i] = sqrt(self.__anorm * (sig2_1 + sig2_2 + sig2_3 + sig2_4))

        return km, sg

    def rodm(self, z):
        """Return the dark matter density

        Keyword arguments:
            z -- redshift
        """
        z1 = 1.0 + z
        ascale = 1.0 / z1
        ascale2 = ascale ** 2.0
        ascale3 = ascale ** 3.0
        return self.__rodm0 / ascale3, 3.0 * self.__rodm0 / ascale2

    def robr(self, z):
        """Return the barionic density.

        Keyword arguments:
            z -- redshift
        """
        z1 = 1.0 + z
        ascale = 1.0 / z1
        ascale3 = ascale ** 3.0
        return self.__robr0 / ascale3

    def roc(self, z):

        return (self.__roc0 / self.H(0) ** 2.0) * self.H(z) ** 2.0

    def age(self, z):
        """Return the age of the Universe for some redshift.

        Keyword arguments:
            z -- redshift
        """
        if(self.__cosmolibImportStatus is True):
            return  self.__lcdmlib.age(z)
        else:
            z1 = 1.0 + z
            ascale = 1.0 / z1
            ascale3 = ascale ** 3.0
            fct = self.__omegalm * ascale3
            return 6.522916e+09 * log(sqrt(fct) + sqrt(fct + 1.0)) / self.__hsl

    def omegamz(self, z):
        om = (self.H(z) / self.H(0)) ** 2.0 - self.__omegal

        om = om * (1.0 / (1.0 + z))

        return om

    def setCosmologicalParameter(self, omegam, omegab, omegal, h):
        """Set the cosmological parameters

        """
        self.__omegab = omegab
        self.__omegam = omegam
        self.__omegal = omegal
        self.__h = h
        return True

    def getCosmologicalParameter(self):
        """Return the cosmological parameter
        """
        return self.__omegab, self.__omegam, self.__omegal, self.__h

    def getDeltaC(self):
        """Return the critical density
        """
        return self.__deltac

    def getTilt(self):
        return self.__tilt

    def getRobr0(self):
        """Return the barionic matter density at the present day.
        """
        return self.__robr0

    def getRodm0(self):
        """Return the dark matter density at the present day.
        """
        return self.__rodm0
