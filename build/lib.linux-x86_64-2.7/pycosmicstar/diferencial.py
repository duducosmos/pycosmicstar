#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import division, absolute_import

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"
"""

    This file is part of pystar.
    copyright : Eduardo dos Santos Pereira
    31 mar. 2011.

    pystar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.
    pystar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

Fornce as funcoes locate(xx,n,x), dfridr(func,x,h,err) e
int_simples(func,a,b,dx =0.001)
"""

from numpy import zeros
from numpy import float64 as Float64
from numpy import abs


def locate(xx, n, x):
    """Localiza a posicao de dado ponto a partir de dois adjacentes.

argumentos:  func --- funcao ou tabela de entrada
xx   --- tabela de entrada
n    --- numero de pontos da tabela
x    --- valor de x que se deseja determinar y
j    --- posicao de saida
"""
    jl = 0
    ju = n + 1
    while(ju - jl > 1):
        jm = int((ju + jl) / 2)
        if(xx[n] > xx[1] and x > xx[jm]):
            jl = jm
        else:
            ju = jm
    return jl


def dfridr(func, x, h, err):
    '''Fornece a derivada de y em relacao a x.
           argumentos:  func --- funcao a ser integrada
                        x    --- dlog10 m ou z
                        h    --- passo para a diferencicao
                        err  --- parametro interno de erro da function
    '''
    CON = 1.4
    CON2 = CON * CON
    #BIG = 1.0e30
    NTAB = 10
    #SAFE = 2.0
    dfit = 0.0
    a = zeros((NTAB, NTAB), type=Float64)
    if(h == 0):
        print('h tem que ser diferente de zero')
        return
    hh = h
    a[0, 0] = (func(x + hh) - func(x - hh)) / (2.0 * hh)
    for i in range(1, NTAB):
        hh = hh / CON
        a[0, i] = (func(x + hh) - func(x - hh)) / (2.0 * hh)
        fac = CON2
        for j in range(1, i):
            a[j, i] = (a[j - 1, i] * fac - a[j - 1, i - 1]) / (fac - 1.0)
            fac = CON2 * fac
            a1 = a[j, i] - a[j - 1, i]
            a2 = a[j, i] - a[j - 1, i - 1]
            errt = max(abs(a1), abs(a2))
            if (errt >= err):
                err = errt
                dfit = a[j, i]
                return dfit
    return None
