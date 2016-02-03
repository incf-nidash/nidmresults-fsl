#!/usr/bin/env python
"""
Test NIDM FSL export tool installation


@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2015
"""
import unittest
from subprocess import check_call


class TestInstall(unittest.TestCase):

    def test_install(self):
        """
        Test: Check that nidmfsl was installed properly
        """
        check_call("python setup.py install; nidmfsl -h", shell=True)

if __name__ == '__main__':
    unittest.main()
