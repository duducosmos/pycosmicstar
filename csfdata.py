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


"""Cosmic Star Formation Rate Observational Data

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

from numpy import array, loadtxt, arange
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from pystar.lcdmcosmology import lcdmcosmology
from pystar.cosmicstarformation import cosmicstarformation


class csfdata:

    def __init__(self):

        """"Data description: The first collunm is the redshift the
        3th and 4th are the top and up errors in z.
        The 2th collunm is the cosmic star formation rate (csfr) and the
        5th and 6th collunm are the top and up errors in csfr."""

        arq = open("./data/hopkins_2004.dat", 'r')
        self.data = loadtxt(arq, delimiter=",")
        self.tauBestFit = 2.29

    def csfredshift(self):
        """Return the redshift and the CSFR from
        observational data
        """
        return self.data[:, 0], self.data[:, 1]

    def errorData(self):
        """Return the asymetric errors in the redshif and CSFR
        respectively
        """
        xerr = array([self.data[:, 0] - self.data[:, 2],
                self.data[:, 3] - self.data[:, 0]])
        yerr = array([self.data[:, 1] - self.data[:, 4],
                      self.data[:, 5] - self.data[:, 1]])
        return xerr, yerr

    def plotObservationalData(self):
        x, y = self.csfredshift()
        xerr, yerr = self.errorData()
        plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='.')
        plt.yscale('log')
        plt.show()

    def funcMinimize(self, p, x, yn):
        tau = p[0]

        if(tau > 3.5):
            return yn * 1e9
        elif(tau < 1.0):
            return yn * 1e9

        myCSFR = cosmicstarformation(cosmology=lcdmcosmology,
                                     tau=tau)

        yt = array([myCSFR.cosmicStarFormationRate(zi) for zi in x])

        return yn - yt

    def funcMinimizeTinker(self, p, x, yn):
        tau = p[0]

        #DHs = [200, 300, 400, 600, 800, 1200, 1600, 2400]

        myCSFR = cosmicstarformation(cosmology=lcdmcosmology,
                                     tau=tau,
                                     massFunctionType="TK",
                                     delta_halo=400)

        yt = array([myCSFR.cosmicStarFormationRate(zi) for zi in x])

        return yn - yt

    def fitTauPereiraMiranda(self):
        x, yn = self.csfredshift()
        xerr, yerr = self.errorData()

        p = 2.29090099
        p, sucess = leastsq(self.funcMinimize, p, args=(x, yn))
        return p[0]

    def fitTauTinker(self):
        x, yn = self.csfredshift()
        xerr, yerr = self.errorData()

        p = [0.5]
        p, sucess = leastsq(self.funcMinimizeTinker, p, args=(x, yn))
        return p

    def plotTauBestFit(self):

        myCSFR = cosmicstarformation(cosmology=lcdmcosmology,
                                     tau=self.tauBestFit)
        x, yn = self.csfredshift()
        xerr, yerr = self.errorData()

        z = arange(0, 7.1, 0.1)

        yt = array([myCSFR.cosmicStarFormationRate(zi) for zi in z])

        plt.plot(z, yt, label=r"PM - $\tau = $" + str(self.tauBestFit)
                        + " Gyr")
        plt.errorbar(x, yn, yerr=yerr, xerr=xerr, fmt='.')
        plt.ylabel(r"$\dot{\rho}_{*}$(M$_{\odot}$Mpc$^{-3}$yr$^{-1}$)")
        plt.xlabel(r"$z$")
        plt.yscale("log")
        plt.legend(loc=4)

        plt.show()

    def plotCSFR(self):
        myCSFR_ST = cosmicstarformation(cosmology=lcdmcosmology)
        myCSFR_TK = cosmicstarformation(cosmology=lcdmcosmology,
                                        massFunctionType="TK",
                                        tau=0.64014971,
                                        delta_halo=400)
        myCSFR_TK.setDeltaHTinker(2400)
        z = arange(0, 7.1, 0.1)

        csfrST = array([myCSFR_ST.cosmicStarFormationRate(zi) for zi in z])

        csfrTK = array([myCSFR_TK.cosmicStarFormationRate(zi) for zi in z])

        x, yn = self.csfredshift()
        xerr, yerr = self.errorData()

        plt.plot(z, csfrST, label="ST")
        plt.plot(z, csfrTK, label="TK")
        plt.errorbar(x, yn, yerr=yerr, xerr=xerr, fmt='.')
        plt.legend(loc=4)
        plt.show()


if(__name__ == "__main__"):
    myCSF = csfdata()
    myCSF.plotCSFR()
    #print((myCSF.fitTauTinker()))
