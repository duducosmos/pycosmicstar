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


"""Cosmological Dark Halos History
From the formalism of Reed et al (MNRAS, 346, 565-572, 2003)
it is calculated the mass fraction of dark matter halos.
The code obtain the mass density of dark halos and the fraction
of brions into structures as a function of the time.
Here is used the Transfer function from Efstathiou, Bond & White
-- (MNRAS, 258, 1P, 1992).
The current version it is assumed the normalization of WMAP (withou
gravitational waves) adapted from Eisenstein e Hu (ApJ 511, 5, 1999) that
in the way that return  sigma_8 = 0,84.
The fraction of mass of dark halos is obtained by the work of
Sheth e Tormen (MNRAS 308, 119, 1999).
All models consider Omega_Total = Omega_M + Omega_L = 1,0

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

from numpy import sqrt, pi, log, log10, exp, array, abs
import scipy.interpolate as spint
from scipy.interpolate import InterpolatedUnivariateSpline as spline

from .structuresabstract import structuresabstract

import sys
pyversion = sys.version_info
if pyversion[0] >= 3:
    from . import filedict
else:
    print("Importing filedict for python2.7")
    from . import filedict_old as filedict

import os
from .diferencial import dfridr, locate

from .paralleloverlist import parallel_list


class structures(structuresabstract):
    """This class was contructed based in the like Press-Schechter formalism
    that provides characteristis like numerical density of dark matter halos
    into the range m_h, m_h + dm_h, the fraction of barionic matter,
    and,  the accretion rate of barions into structures and the total number
    of dark halos.

    The models used to develop this class was presented for the first time
    in the article of Pereira and Miranda (2010) - (MNRAS, 401, 1924, 2010).

    The cosmologic background model is passed as a instance parameter:
        cosmology

    Keyword arguments:
        lmin -- (default 6.0) log10 of the minal mass of the dark halo
                            where it is possible to have star formation.

        zmax -- (defaul 20.0) - the maximum redshift to be considered

        omegam -- (default 0.24) - The dark matter parameter

        omegab -- (default 0.04) - The barionic parameter

        omegal -- (default 0.73) - The dark energy parameter

        h -- (default 0.7) - The h of the Hubble constant (H = h * 100)
    """

    def __init__(self, cosmology, lmin=6.0, zmax=20.0,
                omegam=0.24, omegab=0.04, omegal=0.73, h=0.7,
                cacheDir=None, cacheFile=None, massFunctionType="ST",
                delta_halo=200):

        self._cosmology = cosmology(omegam, omegab, omegal, h)

        if(cacheDir is None):
            self._cacheDir = self. __creatCachDiretory()[0]
        else:
            self._cacheDir = cacheDir

        if(cacheFile is None):
            if(massFunctionType == "TK"):

                cacheFile = str(self._cacheDir) + "/structures_cache_"\
                  + massFunctionType + str(delta_halo) + "_" + "_" +\
                   str(omegab) + "_" \
                  + str(omegam) + "_" +\
                   str(omegal) + "_ " \
                  + str(h) + "_" + str(lmin) + "_" + str(zmax)

            else:
                cacheFile = str(self._cacheDir) + "/structures_cache_"\
                      + massFunctionType + "_" + str(omegab) + "_" \
                      + str(omegam) + "_" + str(omegal) + "_ " \
                      + str(h) + "_" + str(lmin) + "_" + str(zmax)
        else:
            cacheFile = str(self._cacheDir) + cacheFile

        self._cache_dict = filedict.FileDict(filename=cacheFile + ".cache")

        self.__mmin = 1.0e+4
        self.__mmax = 1.0e+18
        self.__lmax = log10(self.__mmax / 10.0)
        self._zmax = zmax
        self.__lmin = lmin
        self.__deltac = self._cosmology.getDeltaC()
        self.__pst = 0.3

        h2 = h * h
        h2om = h2 * omegam
        #h2br = h2 * omegab
        self.__ut = 1.0 / 3.0
        self.__nr = 14000
        self.__ct0 = 4.0 * pi
        self.__ct1 = self.__ct0 * 2.76e+11 / 3.0
        self.__ct2 = self.__ct1 * h2om
        self.__ast1 = 0.322
        self.__ast2 = 0.707
        self.__pst = 0.3
        self.__tilt2 = self._cosmology.getTilt() / log(10.0)
        self.__ctst = self.__ast1 * sqrt(2.0 * self.__ast2 / pi)

        self.__massFunctionType = massFunctionType
        self.__delta_halo = delta_halo

        self.__massFunctionDict = {"ST": self.__massFunctionST,
                                   "TK": self.__massFunctionTinker,
                                   "PS": self.__massFunctionPressSchechter,
                                   "JK": self.__massFunctionJenkins,
                                   "W": self.__massFunctionW
                                   }

        self.__startingSigmaAccretion()

    def __creatCachDiretory(self):
        HOME = os.path.expanduser('~')
        if not os.path.exists(HOME + '/.cosmicstarformation'):
            print(('Creating .cosmicstarformation cache diretory in %s'
                     % HOME))
            os.makedirs(HOME + '/.cosmicstarformation')
        return os.path.expanduser('~') + '/.cosmicstarformation', True

    def __ifSigmaNotInCache(self):
        """Calculate the values necessaries to initialize the
        numerical function of sigma
        """
        numk = 1000.0
        kscale = self.__mmax / self.__mmin
        kls = log10(kscale)
        numk = numk * kls
        kls1 = kls / numk
        deltaz = self._zmax / (numk)

        def CalculaKm(i):
            kmass = (10.0 ** ((i + 1) * kls1)) * self.__mmin
            return kmass

        def CalculaScale(i):
            scale = (CalculaKm(i) / self.__ct2) ** self.__ut
            return scale

        self.__kmass = array([CalculaKm(i) for i in range(int(numk))])
        self.__scale = array([CalculaScale(i) for i in range(int(numk))])
        self.__zred = array([self._zmax - i * deltaz
                                  for i in range(int(numk))])

        e, f = self._cosmology.sigma(self.__kmass)

        self.__km = array([ei for ei in e])
        self.__sg = array([FI for FI in f])

        self.__t_z = parallel_list(self._cosmology.age, self.__zred)

        self.__d_c2 = parallel_list(self.__deltaCz, self.__zred)

        self.__rdm2 = parallel_list(self.__rodmz, self.__zred)

        self.__rbr2 = parallel_list(self._cosmology.robr, self.__zred)

        #self.__t_z = array([self._cosmology.age(zi) for zi in self.__zred])

        #self.__d_c2 = array([
            #self.__deltac / self._cosmology.growthFunction(zi)
                            #for zi in self.__zred
                            #])
        #self.__rdm2 = array([self._cosmology.rodm(zi)[0]
                            #for zi in self.__zred
                            #])
        #self.__rbr2 = array([self._cosmology.robr(zi)
                                #for zi in self.__zred])

    def __rodmz(self, z):
        return self._cosmology.rodm(z)[0]

    def __deltaCz(self, z):
        return self.__deltac / self._cosmology.growthFunction(z)

    def __startingSigmaAccretion(self):
        """
        Verify if the values are in cache
        """

        try:

            self.__km = self._cache_dict['km']
            self.__scale = self._cache_dict['scale']
            self.__zred = self._cache_dict['zred']
            self.__sg = self._cache_dict['sg']
            self.__t_z = self._cache_dict['t_z']
            self.__d_c2 = self._cache_dict['d_c2']
            self.__rdm2 = self._cache_dict['rdm2']
            self.__rbr2 = self._cache_dict['rbr2']
            self._abt2 = self._cache_dict['abt2']
            self._ascale = self._cache_dict['ascale']
            self._tck_ab = self._cache_dict['tck_ab']

            print("Data in Cache")

        except:
            self.__ifSigmaNotInCache()
            self.__startBarionicAccretionRate()
            self.__cachingAtribut()

    def __cachingAtribut(self):
        """Caching the values
        """
        self._cache_dict['km'] = self.__km
        self._cache_dict['scale'] = self.__scale
        self._cache_dict['zred'] = self.__zred
        self._cache_dict['sg'] = self.__sg
        self._cache_dict['t_z'] = self.__t_z
        self._cache_dict['d_c2'] = self.__d_c2
        self._cache_dict['rdm2'] = self.__rdm2
        self._cache_dict['rbr2'] = self.__rbr2
        self._cache_dict['abt2'] = self._abt2
        self._cache_dict['ascale'] = self._ascale
        self._cache_dict['tck_ab'] = self._tck_ab

    def massFunction(self, lm, z):
        """Return the mass function of dark halos.

        Keyword arguments:
            lm -- log10 of the mass of the dark halo
            z -- redshift
        """
        try:
            return self.__massFunctionDict[self.__massFunctionType](lm, z)
        except:
            raise NameError("No Defined Mass Function")

    def __massFunctionJenkins(self, lm, z):
        """Return the mass function of Jenkins et al. (2003).
         Keyword arguments:
            lm -- log10 of the mass of the dark halo
            z -- redshift
        """

        gte = self._cosmology.growthFunction(z)
        rdmt, drdmt = self._cosmology.rodm(z)
        step = lm / 2.0e+1
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        sgm = self.fstm(lm)
        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)
        sigma1 = self.__deltac / (sgm * gte)
        #sigma2 = sigma1 ** 2.0

        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)

        fst = 0.315 * exp(- abs(log(sigma1) + 0.61) ** 3.8)

        frst = (rdmt / kmass ** 2.0) * fst * abs(dsgm_dlgm) / sgm
        dn_dm = frst
        return dn_dm

    def __massFunctionPressSchechter(self, lm, z):
        """Return the value of Press-Schechter (1974) mass function.
        Keyword arguments:
            lm -- log10 of the mass of the dark halo
            z -- redshift
        """

        gte = self._cosmology.growthFunction(z)
        rdmt, drdmt = self._cosmology.rodm(z)
        step = lm / 2.0e+1
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        sgm = self.fstm(lm)
        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)
        sigma1 = self.__deltac / (sgm * gte)
        sigma2 = sigma1 ** 2.0
        fst = sqrt(2.0 / pi) * (sigma1) * exp(-0.5 * sigma2)
        frst = (rdmt / kmass ** 2.0) * fst * abs(dsgm_dlgm) / sgm
        dn_dm = frst
        return dn_dm

    def __massFunctionTinker(self, lm, z):
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

        A_0 = A_func(self.__delta_halo)
        a_0 = a_func(self.__delta_halo)
        b_0 = b_func(self.__delta_halo)
        c_0 = c_func(self.__delta_halo)

        A = A_0 * (1 + z) ** (-0.14)
        a = a_0 * (1 + z) ** (-0.06)
        alpha = exp(-(0.75 / log(self.__delta_halo / 75)) ** 1.2)
        b = b_0 * (1 + z) ** (-alpha)
        c = c_0

        gte = self._cosmology.growthFunction(z)
        rdmt, drdmt = self._cosmology.rodm(z)
        step = lm / 2.0e+1
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        sgm = self.fstm(lm)
        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)
        sigma1 = self.__deltac / (sgm * gte)
        sigma2 = sigma1 ** 2.0

        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)

        fst = A * ((sgm / b) ** (-a) + 1) * exp(-c / sgm ** 2.0)

        frst = (rdmt / kmass ** 2.0) * fst * abs(dsgm_dlgm) / sgm
        dn_dm = frst
        return dn_dm

    def __massFunctionW(self, lm, z):
        # LANL fitting function - Warren et al. 2005, astro-ph/0506395, eqtn. 5
        A = 0.7234
        a = 1.625
        b = 0.2538
        c = 1.1982

        gte = self._cosmology.growthFunction(z)
        rdmt, drdmt = self._cosmology.rodm(z)
        step = lm / 2.0e+1
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        sgm = self.fstm(lm)
        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)
        sigma1 = self.__deltac / (sgm * gte)
        sigma2 = sigma1 ** 2.0

        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)

        fst = A * ((sigma1 ** (-a)) + b) * exp(-c / sigma2)

        frst = (rdmt / kmass ** 2.0) * fst * abs(dsgm_dlgm) / sgm
        dn_dm = frst
        return dn_dm

    def __massFunctionST(self, lm, z):
        """Return the mass function of dark halos of
        Sheth e Tormen (MNRAS 308, 119, 1999).

        Keyword arguments:
            lm -- log10 of the mass of the dark halo
            z -- redshift
        """
        gte = self._cosmology.growthFunction(z)
        #gte2 = self._cosmology.dgrowth_dt(z)
        rdmt, drdmt = self._cosmology.rodm(z)
        step = lm / 2.0e+1
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        sgm = self.fstm(lm)
        dsgm_dlgm = dfridr(self.fstm, lm, step, err=0.0)
        sigma1 = self.__deltac / (sgm * gte)
        sigma2 = sigma1 ** 2.0
        expn = exp(-self.__ast2 * sigma2 / 2.0)
        fst = self.__ctst * sigma1 * \
            (1.0 + (1.0 / (sigma2 * self.__ast2)) ** self.__pst) * expn
        frst = (rdmt / kmass ** 2.0) * fst * abs(dsgm_dlgm) / sgm
        dn_dm = frst
        return dn_dm

    def fstm(self, lm):
        '''Numerical function that return the value of sigm that
        will be used by dfridr to calculate d_sigma_dlog10(m).

        Keyword arguments:
            lm -- log10 of the mass of dark halo
        '''
        j = locate(self.__km, len(self.__km) - 1, lm)
        return self.__sg[j]

    def __fmassM(self, lm, z):
        """Return the mass function of dark halos multiplied by Mass -
        Sheth e Tormen (MNRAS 308, 119, 1999).
        """
        kmsgm = lm
        kmass = 10.0 ** (kmsgm)
        frst = self.massFunction(lm, z) * kmass
        kmassa2 = self.__tilt2 * kmass
        mdn_dm = kmassa2 * frst
        return mdn_dm

    def halos_n(self, z):
        """Return the integral of the mass function of dark halos multiplied
        by mass in the range of log(M_min) a log(M_max)

        Keyword arguments:
            z -- redshift
        """

        fmassM = lambda lm: self.__fmassM(lm, z)

        deltal = (self.__lmax - self.__lmin) / 50.0

        Lm = [self.__lmin + i * deltal for i in range(50)]

        Fm = [fmassM(lm) for lm in Lm]

        tck = spint.splrep(Lm, Fm)
        Inte = spint.splint(self.__lmin, self.__lmax, tck)
        return Inte

    def fbstruc(self, z):
        """Return the faction of barions into structures

        Keyword arguments:
            z -- redshift
        """
        rdm, drdm_dt = self._cosmology.rodm(z)
        fb = self.halos_n(z) / rdm
        return fb

    def numerical_density_halos(self, z):
        """Return the numerial density of dark halos
        within the comove volume

        Keyword arguments:
            z- redshift
        """

        deltal = (self.__lmax - self.__lmin) / 50.0
        Lm = [self.__lmin + i * deltal for i in range(50)]

        Fm = [self.massFunction(lm, z) for lm in Lm]

        tck = spint.splrep(Lm, Fm)
        Inte = spint.splint(self.__lmin, self.__lmax, tck)
        return Inte

    def abt(self, a):
        """Return the accretion rate of barionic matter, as
        a function of scala factor, into strutures.

        Keyword arguments:
            a -- scala factor (1.0 / (1.0 + z))
        """
        i = locate(self._ascale, len(self._ascale) - 1, a)
        return self._abt2[i]

    def __startBarionicAccretionRate(self):

        np = 1000
        deltaz = self._zmax / float(np)

        z = [self._zmax - i * deltaz for i in range(np)]
        z.append(0)
        z = array(z)
        fbt2 = array([self.fbstruc(zi) for zi in z])
        ascale = array([1.0 / (1.0 + zi) for zi in z])
        self._ascale = ascale

        tck = spint.splrep(ascale, fbt2)
        ab3 = spint.splev(ascale, tck, der=1)

        def a5(z, i):
            a = 1.0 / (1.0 + z)
            a2 = a * a
            a3 = -1.0 * ab3[i] * a2
            a4 = a3
            a5 = self._cosmology.getRobr0() * abs(a4) \
                 / self._cosmology.dt_dz(z)
            return a5

        self._abt2 = array([a5(z[i], i) for i in range(z.size)])
        self._tck_ab = spint.splrep(self._ascale, self._abt2)

    def getCacheDir(self):
        """Return True and cache name if the cache directory existe
        and false else.
        """
        if(self._cacheDir is not None):
            return True, self._cacheDir
        else:
            return False

    def setDeltaHTinker(self, delta_halo):
        if(self.__massFunctionType == "TK"):
            self.__delta_halo = delta_halo
            return True
        else:
            return False

    def getDeltaHTinker(self):
        if(self.__massFunctionType == "TK"):
            return self.__delta_halo
        else:
            raise NameError("TinkerNotDefined")
