#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

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
"""

from numpy import sqrt, pi, log, log10, exp, array, abs
import scipy.interpolate as spint
from structuresPS import structuresPS

import filedict
import os
from diferencial import dfridr, locate


class structures(structuresPS):
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
                cacheDir=None):

        self.__cosmology = cosmology(omegam, omegab, omegal, h)

        if(cacheDir is None):
            self.__cacheDir = self. __creatCachDiretory()[0]
        else:
            self.__cacheDir = cacheDir

        arq = str(self.__cacheDir) + "/structures_cache_" + str(omegab) + "_" \
              + str(omegam) + "_" + str(omegal) + "_" \
              + str(h) + "_" + str(lmin) + "_" + str(zmax)

        self.__cache_dict = filedict.FileDict(filename=arq + ".cache")

        self.__mmin = 1.0e+4
        self.__mmax = 1.0e+18
        self.__lmax = log10(self.__mmax / 10.0)
        self.__zmax = zmax
        self.__lmin = lmin
        self.__deltac = self.__cosmology.getDeltaC()
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
        self.__tilt2 = self.__cosmology.getTilt() / log(10.0)
        self.__ctst = self.__ast1 * sqrt(2.0 * self.__ast2 / pi)

        self.__startingSigmaAccretion()

    def __creatCachDiretory(self):
        HOME = os.path.expanduser('~')
        if not os.path.exists(HOME + '/.scc-strarsFormation'):
            print(('Creating .scc-strarsFormation cache diretory in %s' % HOME))
            os.makedirs(HOME + '/.scc-strarsFormation')
        return os.path.expanduser('~') + '/.scc-strarsFormation', True

    def __ifSigmaNotInCache(self):
        """Calculate the values necessaries to initialize the
        numerical function of sigma
        """
        numk = 1000.0
        kscale = self.__mmax / self.__mmin
        kls = log10(kscale)
        numk = numk * kls
        kls1 = kls / numk
        deltaz = self.__zmax / (numk)

        def CalculaKm(i):
            kmass = (10.0 ** ((i + 1) * kls1)) * self.__mmin
            return kmass

        def CalculaScale(i):
            scale = (CalculaKm(i) / self.__ct2) ** self.__ut
            return scale

        self.__kmass = array([CalculaKm(i) for i in range(int(numk))])
        self.__scale = array([CalculaScale(i) for i in range(int(numk))])
        self.__zred = array([self.__zmax - i * deltaz
                                  for i in range(int(numk))])

        e, f = self.__cosmology.sigma(self.__kmass)

        self.__km = array([ei for ei in e])
        self.__sg = array([FI for FI in f])
        self.__t_z = array([self.__cosmology.age(zi) for zi in self.__zred])
        self.__d_c2 = array([
            self.__deltac / self.__cosmology.growthFunction(zi)
                            for zi in self.__zred
                            ])
        self.__rdm2 = array([self.__cosmology.rodm(zi)[0]
                            for zi in self.__zred
                            ])
        self.__rbr2 = array([self.__cosmology.robr(zi) for zi in self.__zred])

    def __startingSigmaAccretion(self):
        """
        Verify if the values are in cache
        """

        try:

            self.__km = self.__cache_dict['km']
            self.__scale = self.__cache_dict['scale']
            self.__zred = self.__cache_dict['zred']
            self.__sg = self.__cache_dict['sg']
            self.__t_z = self.__cache_dict['t_z']
            self.__d_c2 = self.__cache_dict['d_c2']
            self.__rdm2 = self.__cache_dict['rdm2']
            self.__rbr2 = self.__cache_dict['rbr2']
            self.__abt2 = self.__cache_dict['abt2']
            self.__ascale = self.__cache_dict['ascale']
            self.__tck_ab = self.__cache_dict['tck_ab']

            print("Data in Cache")

        except:
            self.__ifSigmaNotInCache()
            self.__startBarionicAccretionRate()
            self.__cachingAtribut()

    def __cachingAtribut(self):
        """Caching the values
        """
        self.__cache_dict['km'] = self.__km
        self.__cache_dict['scale'] = self.__scale
        self.__cache_dict['zred'] = self.__zred
        self.__cache_dict['sg'] = self.__sg
        self.__cache_dict['t_z'] = self.__t_z
        self.__cache_dict['d_c2'] = self.__d_c2
        self.__cache_dict['rdm2'] = self.__rdm2
        self.__cache_dict['rbr2'] = self.__rbr2
        self.__cache_dict['abt2'] = self.__abt2
        self.__cache_dict['ascale'] = self.__ascale
        self.__cache_dict['tck_ab'] = self.__tck_ab

    def funcMassST(self, lm, z):
        """Return the mass function of dark halos of
        Sheth e Tormen (MNRAS 308, 119, 1999).

        Keyword arguments:
            lm -- log10 of the mass of the dark halo
            z -- redshift
        """
        gte = self.__cosmology.growthFunction(z)
        #gte2 = self.__cosmology.dgrowth_dt(z)
        rdmt, drdmt = self.__cosmology.rodm(z)
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
        frst = self.funcMassST(lm, z) * kmass
        kmassa2 = self.__tilt2 * kmass
        mdn_dm = kmassa2 * frst
        return mdn_dm

    def halos_n(self, z):
        """Return the integral of the mass function of dark halos multiplied
        by mass in the range of log(M_min) a log(M_max)

        Keyword arguments:
            z -- redshift
        """

        fmassM = lambda ln: self.__fmassM(lm, z)

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
        rdm, drdm_dt = self.__cosmology.rodm(z)
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

        Fm = [self.funcMassST(lm, z) for lm in Lm]

        tck = spint.splrep(Lm, Fm)
        Inte = spint.splint(self.__lmin, self.__lmax, tck)
        return Inte

    def abt(self, a):
        """Return the accretion rate of barionic matter, as
        a function of scala factor, into strutures.

        Keyword arguments:
            a -- scala factor (1.0 / (1.0 + z))
        """
        i = locate(self.__ascale, len(self.__ascale) - 1, a)
        return self.__abt2[i]

    def __startBarionicAccretionRate(self):

        np = 1000
        deltaz = self.__zmax / float(np)

        z = [self.__zmax - i * deltaz for i in range(np)]
        z.append(0)
        z = array(z)
        fbt2 = array([self.fbstruc(zi) for zi in z])
        ascale = array([1.0 / (1.0 + zi) for zi in z])
        self.__ascale = ascale

        tck = spint.splrep(ascale, fbt2)
        ab3 = spint.splev(ascale, tck, der=1)

        def a5(z, i):
            a = 1.0 / (1.0 + z)
            a2 = a * a
            a3 = -1.0 * ab3[i] * a2
            a4 = a3
            a5 = self.__cosmology.getRobr0() * abs(a4) \
                 / self.__cosmology.dt_dz(z)
            return a5

        self.__abt2 = array([a5(z[i], i) for i in range(z.size)])
        self.__tck_ab = spint.splrep(self.__ascale, self.__abt2)

    def getCacheDir(self):
        """Return True if the cache directory existe and false else.
        """
        if(self.__cacheDir is not None):
            return True
        else:
            return False
