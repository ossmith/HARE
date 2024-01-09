#!/usr/in/env python
# sigtest_test.py

# Example command: python sigtest_test.py

from pandas import read_table
import numpy as np
import unittest
import sigtest
import os
from sys import argv
from sys import exit
import glob
# from cv2 import imread
# import filecmp

working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestSigtest(unittest.TestCase):

    def setUp(self):
        self.passIntersections = f"{working_dir}/input/sigtest_pass.intersections"
        data = read_table(self.passIntersections, sep="\t", header=0)
        self.simulations = data[data["category"] == "simulation"]
        self.test_value = data[data["category"] == "test_set"]
        self.p_values = [0.2231, 0.3031, 0.0013, 0.0062, 0.0019, 0.8240, 0.5240, 0.9458, 1.0, 0.6600, 0.0000]

    def tearDown(self):
        # Clean up directory and clear cache
        for f in glob.glob(f"{working_dir}/test.*"):
            os.remove(f)
        # filecmp.clear_cache()

    def test_sigtest(self):
        # print("test_sigtest")
        passStats = [500, 100, 1.0298275760674873e-06, 1.4043103310011189e-06, 0.08]
        (sdf_ignore, tdf_ignore, testStats) = sigtest.getStats(self.passIntersections)
        self.assertEqual(testStats[1:], passStats)

    def test_weibull(self):
        # print("test_weibull")
        weibullResults = sigtest.weibullTest(self.simulations, self.test_value, 74899)
        weibullPass = [5.0987075897638645, 1.117635750036753e-06, 0.04062655117962466] # c, x, weibull_p
        self.assertEqual(weibullPass, weibullResults)

    def test_plot(self):
        # print("test_plot")
        plot_pass = f"{working_dir}/input/plot_pass.png"
        out_stem = f"{working_dir}/test"
        simDF = read_table(f"{working_dir}/input/plot_inS.txt", header=None, names=["int_per_bp"], dtype="float")
        testDF = read_table(f"{working_dir}/input/plot_inT.txt", header=0)
        checkRun = sigtest.plot(simDF, testDF, "sigtest_pass", 0, out_stem)
        self.assertTrue(checkRun) # Basic test to check if it ran, does not currently check that files are the same

        # Check equivalent files with cv2
        # plotP = imread(plot_pass)
        # plotT = imread(f"{out_stem}.png")
        # self.assertTrue(np.all(plotP==plotT))

        # Check equivalent files with filecmp
        # self.assertTrue(filecmp.cmp(f"{out_stem}.png", plot_pass, shallow=True)) # Use filecmp; only works if matplotlib versions are the same

    def test_bhCorrection(self):
        pass_bhp = [0.0000, 0.00715, 0.006966666666666666, 0.01705, 0.49082, 0.5556833333333333, 0.8234285714285715, 0.9075000000000001, 1.00, 1.00, 1.00]
        out_bhp = sigtest.bhCorrection(self.p_values)
        self.assertTrue(out_bhp == pass_bhp)

if __name__ == '__main__':
    unittest.main()
