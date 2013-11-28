#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

"""
Abstract Class of cosmological models.
"""


class Cosmology:

    def dt_dz(self, z):
        "Return the relation between the cosmic time and the redshift"
        return NotImplementedError('I need to be implemented!')

    def dr_dz(self, z):
        "Return the comove-radii"
        return NotImplementedError('I need to be implemented!')

    def dV_dz(self, z):
        "Return the comove volume"
        return NotImplementedError('I need to be implemented!')

    def rhoDarkMatter(z):
        "Return the Dark Matter Density"
        return NotImplementedError('I need to be implemented!')

    def rhoBarionicMatter(z):
        "Return the Barionic Matter Density"
        return NotImplementedError('I need to be implemented!')

    def H(z):
        "Return the Hubble Parameter"
        return NotImplementedError('I need to be implemented!')


