#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

"""
Abstract Class of cosmological models.
"""


class cosmology:

    def dt_dz(self, z):
        """Return the relation between the cosmic time and the redshift
        """
        raise NotImplementedError('I need to be implemented!')

    def dr_dz(self, z):
        """Return the comove-radii"""
        raise NotImplementedError('I need to be implemented!')

    def dV_dz(self, z):
        """Return the comove volume"""
        raise NotImplementedError('I need to be implemented!')

    def rodm(self, z):
        """Return the Dark Matter Density"""
        raise NotImplementedError('I need to be implemented!')

    def robr(self, z):
        """Return the Barionic Matter Density"""
        raise NotImplementedError('I need to be implemented!')

    def H(self, z):
        """Return the Hubble Parameter"""
        raise NotImplementedError('I need to be implemented!')

    def dgrowth_dt(self, z):
        """Return the derivative of growth function of the
         primordial perturbations"""
        raise NotImplementedError('I need to be implemented!')

    def growthFunction(self, z):
        """Return the growth function of the primordial perturbations"""
        raise NotImplementedError('I need to be implemented!')

    #def dsigma2_dk(self, kl):
        #""""Return the integrating of sigma(M,z) for a top-hat filtering.
        #In z = 0 return sigma_8, for z > 0 return sigma(M,z)
        #"""
        #raise NotImplementedError('I need to be implemented!')

    def sigma(self):
        """Return  the variance of the linear density field.
        As pointed out by Jenkis et al. (2001),
        this definition of the mass function has the advantage that it does
        not explicitly depend on redshift, power spectrum or cosmology.
        """
        raise NotImplementedError('I need to be implemented!')

    def age(self, z):
        """Return the age of the Universe for a given z
        """
        raise NotImplementedError('I need to be implemented!')

    def setCosmologicalParameter(self):
        return False

    def getCosmologicalParameter(self):
        return False

    def getTilt(self):
        return False

    def getRobr0(self):
        return False


