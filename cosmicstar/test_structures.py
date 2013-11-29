#!/usr/bin/env python
# *-* Coding: UTF-8 *-*

import unittest

from structures import structures
from lcdmcosmology import lcdmcosmology


class test_structures(unittest.TestCase):

    myStructures = structures(lcdmcosmology)

    def test_funcMassST(self):
        self.assertEquals(self.myStructures.funcMassST(9.0, 1.0),
                          8.4508147954749659e-09)

    def test_fstm(self):
        self.assertEquals(self.myStructures.fstm(6.0), 515.85483654505242)

    def test_halos_n(self):
        self.assertEquals(self.myStructures.halos_n(0.0), 23581044522.84875)

    def test_fbstruc(self):
        self.assertEquals(self.myStructures.fbstruc(0.0), 0.7265184278211807)

    def test_numerical_density_halos(self):
        self.assertEquals(self.myStructures.numerical_density_halos(0.0),
                            6.757090440025199e-05)

    def test_abt(self):
        self.assertEquals(self.myStructures.abt(1.0), 0.0078707904676166528)

    def test_creatCachDiretory(self):
        self.assertTrue(self.myStructures.getCacheDir(),
        "The directory not Exist")



if(__name__ == "__main__"):
    unittest.main()