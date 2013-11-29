#!/usr/bin/env python
# *-* Coding : UTF-8 *-*
import unittest

if(__name__ == "__main__"):
    testsuite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=1).run(testsuite)

