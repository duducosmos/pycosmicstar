#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

"""
Abstract Class of like Press-Schechter formalis
"""


class structuresPS:

    def funcMassST(self, lm, z):
        """Return the mass function of dark halos of
        Sheth e Tormen (MNRAS 308, 119, 1999).
        """
        raise NotImplementedError('I need to be implemented!')

    def fstm(self, lm):
        '''Numerical function that return the value of sigm that
        will be used by dfridr to calculate d_sigma_dlog10(m).'''
        raise NotImplementedError('I need to be implemented!')

    def halos_n(self, z):
        """Return the integral of the mass function of dark halos multiplied
        by mass in the range of log(M_min) a log(M_max)
        """
        raise NotImplementedError('I need to be implemented!')

    def fbstruc(self, z):
        """Return the faction of barions into structures
        """
        raise NotImplementedError('I need to be implemented!')

    def numerical_density_halos(self, z):
        """Return the numerial density of dark halos
        within the comove volume
        """
        raise NotImplementedError('I need to be implemented!')

    def abt(self, a):
        """Return the accretion rate of barionic matter into strutures
        """
        raise NotImplementedError('I need to be implemented!')

    def getCacheDir(self):
        """Return True if the cache directory existe and false else.
        """
        raise NotImplementedError('I need to be implemented!')
