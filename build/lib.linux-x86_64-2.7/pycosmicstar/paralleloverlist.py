#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, absolute_import

__author__ = "Eduardo dos Santos Pereira"
__email__ = "pereira.somoza@gmail.com"
__credits__ = ["Eduardo dos Santos Pereira"]
__license__ = "GPLV3"
__version__ = "1.0.1"
__maintainer__ = "Eduardo dos Santos Pereira"
__status__ = "Stable"
__date__ = '09/12/2013'

"""
"""

import multiprocessing as mpg
from numpy import array

##@file paralleloverlist.py
##@author  Eduardo dos Santos Pereira <pereira.somoza@gmail.com>
##@version 1.1
##@section LICENSE
# This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License as
#published by the Free Software Foundation; either version 2 of
#the License, or (at your option) any later version.
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details at
#http://www.gnu.org/copyleft/gpl.html
##@section DESCRIPTION


class paralleloverlist:
    """
    ppvector: Parallel Processing Vector
    This program is used to calculate, in parallel, by python module
    multiprocessing, points in vector.
    """

    def __init__(self, func, inputArray):
        '''
        Dmatiz: The dimension of the vector
        func:  function that will run in parallel
        '''
        self.__inputArray = mpg.Array('d', inputArray)
        self.__sizeArray = len(inputArray)
        self.__func = func
        self.__output = mpg.Array('d', inputArray)
        self.__runProcess()

    def getResult(self):
        return array(self.__output)

    def __call__(self):
        return array(self.__output)

    def __Calcula(self, func, k, E, n):
        x = self.__inputArray[k: E + k]
        self.__output[k: E + k] = [func(xi) for xi in x]

    def __acaoParalera(self, n, q, Dmatriz, func, n_process):
        E = Dmatriz // n_process
        k = n * E
        q.put(self.__Calcula(func, k, E, n))

    def __runProcess(self):

        n_process = mpg.cpu_count()
        subprocess = []

        for i in range(n_process):
            q = mpg.Queue()
            p = mpg.Process(target=self.__acaoParalera,
                args=(i, q, self.__sizeArray, self.__func, n_process))
            p.start()
            subprocess.append(p)

        while subprocess:
            subprocess.pop().join()


def parallel_list(func, x):
    result = paralleloverlist(func, x)
    return result()


if(__name__ == "__main__"):
    import time
    import matplotlib.pyplot as plt

    tP = []
    tS = []
    func = lambda x: x ** 2.0
    for i in range(1, 500):
        x = array(list(range(0, i * 100)))
        t1 = time.time()
        parallel_list(func, x)
        t2 = time.time()
        tP.append([i * 100, t2 - t1])
        t3 = time.time()
        [func(xi) for xi in x]
        t4 = time.time()
        tS.append([i * 100, t4 - t3])
    tP = array(tP)
    tS = array(tS)
    plt.plot(tP[:, 1], tP[:, 0], label="Parallel")
    plt.plot(tS[:, 1], tS[:, 0], label="Serial")
    plt.xlabel('time')
    plt.ylabel('size input array')
    plt.legend()
    plt.show()
