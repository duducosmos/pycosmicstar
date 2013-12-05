#!/usr/bin/env python3.3
# *-* Coding: UTF-8 *-*

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"

"""This module contain the Cosmic Star Formation Rate (CSFR) of the work
of:
    Fardal et al. (MNRAS, 379,985,2007) MNRAS, 339,312,2003.
    Hopkins and Beacom, Apj, 651, 142, 2006.
    Volker, Springel, Lars, Hernquist, MNRAS, 339,312,2003.
"""

from numpy import exp


def rho_starF(z, h):
    """Return CSFR by the work of Fardal et al. (MNRAS, 379,985,2007)
    MNRAS, 339,312,2003.
    """
    a = 0.0103
    b = 0.088
    c = 2.4
    d = 2.8
    rho = ((a + b * z) * h) / (1.0 + (z / c) ** d)
    return rho


def rho_starHB(z, h):
    """Return CSFR by the work of Hopkins e Beacom,
    Apj, 651, 142, 2006.
    """
    a = 0.0170
    b = 0.13
    c = 3.3
    d = 5.3
    rho = ((a + b * z) * h) / (1.0 + (z / c) ** d)
    return rho


def rho_starSH(z):
    """Return CSFR by the work of Volker, Springel, Lars and Hernquist
    MNRAS, 339,312,2003.
    """
    #M_{\odot}yr^{-1}Mpc^{-3}
    rhom = 0.15
    alpha = 3.0 / 5.0
    beta = 14.0 / 15.0
    zm = 5.4
    rho = rhom * (beta * exp(alpha * (z - zm))) / \
                (beta - alpha + alpha * exp(beta * (z - zm)))
    return rho
