#!/usr/in/env python
# prerank_test.py

# Example command: python prerank_test.py

import unittest
import prerank
import os
import filecmp
import hareclasses
import argparse
from hare import check_installs
from sys import exit
from sys import argv
import glob
import pandas as pd
from numpy import log10

working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestIntersect(unittest.TestCase):

    def setUp(self):
        self.passParameters = hareclasses.PrerankArgumentContainer(f"{working_dir}/input/gwas_pass.txt", f"{working_dir}/test", "37",
        1000, "all", None, "P", "CHR", "POS", 0.05, 3, 1e-8, 1000, None, None, f"{working_dir}/input/prerank_include.bed",
        "min", False)

        self.xclParameters = hareclasses.PrerankArgumentContainer(f"{working_dir}/input/gwas_pass.txt", f"{working_dir}/test", "37",
        1000, "all", None, "P", "CHR", "POS", 0.05, 3, 1e-8, 1000, None, f"{working_dir}/input/prerank_exclude.bed", None,
        "min", True)

        self.failParameters = hareclasses.PrerankArgumentContainer(f"{working_dir}/input/gwas_pass.txt", f"{working_dir}/test", "37",
        0, "all", None, "P", "CHR", "POS", 1, 3, 1e-8, 1000, None, f"{working_dir}/input/prerank_exclude_fail.bed", f"{working_dir}/input/prerank_include.bed",
        "min", False)

        self.callPeaks = hareclasses.PrerankArgumentContainer(f"{working_dir}/input/gwas_pass.txt", f"{working_dir}/test", "37",
        1000, "all", None, "P", "CHR", "POS", 0.05, 3, 1e-6, 1000, None, None, f"{working_dir}/input/prerank_include.bed",
        "min", True)

        self.biomart_df = pd.read_table(f"{working_dir}/input/prerank_biomart_pass.txt", sep="\t", header=0)

    def tearDown(self):
        # Clean up directory and clear cache
        for f in glob.glob(f"{working_dir}/test.*"):
            os.remove(f)
        filecmp.clear_cache()

    def test_filterPositive(self):
        filtered_pass = f"{working_dir}/input/include_pass.txt"
        incl_out = f"{self.passParameters.output}.filtered.incl"

        filteredDF = prerank.filterPositions(self.passParameters)
        filteredDF.to_csv(incl_out, sep="\t", header=True, index=False)
        self.assertTrue(filecmp.cmp(filtered_pass, incl_out, shallow=False))

    def test_filterNegative(self):
        filtered_xcl = f"{working_dir}/input/exclude_pass.txt"
        xcl_out = f"{self.xclParameters.output}.filtered.xcl"

        filteredDF = prerank.filterPositions(self.xclParameters)
        filteredDF.to_csv(xcl_out, sep="\t", header=True, index=False)
        self.assertTrue(filecmp.cmp(filtered_xcl, xcl_out, shallow=False))

    def test_filter_fail(self):
        with self.assertRaises(ValueError):
            prerank.filterPositions(self.failParameters)

    def test_biomart(self):
        input_df = pd.read_table(f"{working_dir}/input/prerank_pvalDF.txt", sep="\t", header=0)
        biomart_out = f"{self.passParameters.output}.biomart"
        biomart_pass = f"{working_dir}/input/prerank_biomart_pass.txt"

        prerank.annotate(input_df, self.passParameters)
        self.assertTrue(filecmp.cmp(biomart_out, biomart_pass, shallow=False))

    def test_rank(self):
        score_df = pd.read_table(f"{working_dir}/input/prerank_pvalDF.txt", sep="\t", header=0)
        score_out = f"{self.passParameters.output}.rnk"
        score_pass = f"{working_dir}/input/prerank_rnk_pass.rnk"

        prerank.score_and_rank(score_df, self.biomart_df, self.passParameters)
        self.assertTrue(filecmp.cmp(score_out, score_pass, shallow=False))

    def test_call_peaks(self):
        log_df = pd.read_table(f"{working_dir}/input/prerank_pvalDF.txt", sep="\t", header=0)
        log_df["logP"] = log10(log_df["P"])*-1

        called_pass = pd.read_table(f"{working_dir}/input/prerank_calls_pass.txt", sep="\t", header=0)
        pass_sorted = called_pass[["HGNC_SYMBOL","SIG_COUNT"]].sort_values(by="HGNC_SYMBOL", ignore_index=True)

        called_df = prerank.call_peaks(log_df, self.biomart_df, self.callPeaks)
        called_sorted = called_df[["HGNC_SYMBOL","SIG_COUNT"]].sort_values(by="HGNC_SYMBOL", ignore_index=True)

        self.assertTrue(called_sorted.equals(pass_sorted))

if __name__ == '__main__':
    unittest.main()
