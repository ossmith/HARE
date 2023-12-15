#!/usr/in/env python
# sigtest_test.py

from pandas import read_table
import unittest
import sigtest
import os
import filecmp
from sys import argv
from sys import exit
import glob

working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestSigtest(unittest.TestCase):

    def __init__(self, *arguments):
        super().__init__()

        self.test_sigtest()
        self.test_plot()
        self.clean_dir()

    def runTest(self):
        pass

    def test_sigtest(self):
        print("test_sigtest")
        passFile = f"{working_dir}/input/sigtest_pass.intersections"
        passStats = [500, 100, 1.0298275760674873e-06, 1.4043103310011189e-06, 0.08]
        (sdf_ignore, tdf_ignore, testStats) = sigtest.getStats(passFile)
        self.assertEqual(testStats[1:], passStats)
        return

    def test_plot(self):
        print("test_plot")
        plot_pass = f"{working_dir}/input/plot_pass.png"
        out_stem = f"{working_dir}/test"
        simDF = read_table(f"{working_dir}/input/plot_inS.txt", header=None, names=["int_per_bp"], dtype="float")
        testDF = read_table(f"{working_dir}/input/plot_inT.txt", header=0)
        sigtest.plot(simDF, testDF, "sigtest_pass", 0, out_stem)
        self.assertTrue(filecmp.cmp(f"{out_stem}.png", plot_pass, shallow=False))

        return

    def clean_dir(self):
        for f in glob.glob(f"{working_dir}/test.*"):
            os.remove(f)
        filecmp.clear_cache() # Remove filecmp cache
        return

if __name__ == '__main__':
    unittest.main()
