#!/usr/in/env python
# intersect_test.py

# Example command: python intersect_test.py --cache_dir [VEP_CACHE_PATH] --cache_ver [VEP_VERSION]

import unittest
import intersect
import os
import filecmp
import hareclasses
import argparse
from hare import check_installs
from sys import exit
from sys import argv
import glob

parser = argparse.ArgumentParser(description='Automated testing of HARE pipeline.') # Must allow specification of VEP cache
parser.add_argument('--cache_dir', '-c', type=str, help="Specify the VEP cache directory to use. Default is \"$HOME/.vep/\"", required=False, default="$HOME/.vep/")
parser.add_argument('--cache_ver', '-v', type=str, help="Specify the VEP cache version to use. Default is 105. It is recommended by Ensembl that you use the same cache and VEP version.", required=False, default="105")
args = parser.parse_args()

cache_pri = args.cache_dir
cache_dir = os.path.normpath(os.path.expanduser(cache_pri))
# Check directory for VEP cache
if os.path.exists(cache_dir) == False:
    raise FileNotFoundError(f'{cache_dir} does not exist or could not be opened.')
cache_ver = args.cache_ver
working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestIntersect(unittest.TestCase):

    def setUp(self):
        super(TestIntersect, self).setUp()
        self.ref = f"{working_dir}/input/ref_pass.bed"
        self.tools = ["vep", "bedtools"]

        # Pass
        self.passSettings = hareclasses.SettingsContainer(None, None, None, None, True)
        self.passParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_pass.txt", 1e-6, "P", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", self.ref, "37", 3, self.passSettings)

        # Neale
        self.nealeSettings = hareclasses.SettingsContainer(None, None, True, None, None)
        self.nealeParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_neale.txt", 1e-6, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", self.ref, "37", 3, self.nealeSettings)

        # Bolt
        self.boltSettings = hareclasses.SettingsContainer(None, None, None, True, None)
        self.boltParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_bolt.txt", 1e-6, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", self.ref, "37", 3, self.boltSettings)

        # Fail
        self.failParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_fail.txt", 1e-100, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_fail.txt", self.ref, "38", 3, self.passSettings)

        self.settings = [self.passSettings, self.nealeSettings, self.boltSettings]
        self.params = [self.passParameters, self.nealeParameters, self.boltParameters]

        self.snpFiles = [f"{working_dir}/input/snp_pass.snps", f"{working_dir}/input/snp_neale.snps", f"{working_dir}/input/snp_bolt.snps"]

    def tearDown(self):
        # Clean up directory and clear cache
        for f in glob.glob(f"{working_dir}/test.*"):
            os.remove(f)
        filecmp.clear_cache()

    def test_find_files(self):
        intersect.findFiles(self.ref)

    def test_installs(self):
        for t in self.tools:
            check_installs(t)

    # def test_snp2loc(self, testSettings, passParameters):
    #     CHECK SNP MAP
    #     self.assertEqual(intersect.snpToLoc(testSettings, passParameters))
    #
    #     with self.assertRaises(RuntimeError):
    #         intersect.snpToLoc(self) # Mapping failed because mapped_df = 0
    #
    #     return

    def test_header(self):
        defaultHeader = ["CHR", "POS", "P"]
        nealeHeader = ["variant", "pval"]
        boltHeader = ["BP", "SNP", "ALLELE1", "ALLELE0", "P_BOLT_LMM"]
        headersList = [defaultHeader, nealeHeader, boltHeader]
        colList = [["CHR","POS"], ["CHR","POS"], ["CHR","BP"]]

        for h in range(len(headersList)):
            chrResult, posResult = intersect.check_header(self.params[h], self.settings[h], headersList[h], colList[h][0], colList[h][1])
            self.assertEqual(chrResult, colList[h][0])
            self.assertEqual(posResult, colList[h][1])

    def test_header_fail(self):
        with self.assertRaises(KeyError):
            intersect.check_header(self.failParameters, self.passSettings, ["CHR", "POS", "FOO"], "CHR", "POS")

    def test_gwas(self):
        for g in range(len(self.snpFiles)):
            snpLoc = intersect.gwas_import(self.params[g], self.settings[g])
            self.assertTrue(filecmp.cmp(snpLoc, self.snpFiles[g], shallow=False))

    def test_gwas_fail(self):
        with self.assertRaises(RuntimeError):
            intersect.gwas_import(self.failParameters, self.passSettings)

    def test_vep(self):
        snps_out = f"{working_dir}/input/snp_pass.snps"
        # vep_out = f"{argumentClass.output}.features"
        vep_out = f"{self.passParameters.output}.features"
        vep_success = f"{working_dir}/input/vep_pass.features"

        intersect.vep_annotate(snps_out, self.passParameters)
        self.assertTrue(filecmp.cmp(vep_out, vep_success, shallow=False))

    def test_biomart(self):
        vep_out = f"{working_dir}/input/vep_pass.features"
        biomart_out = f"{self.passParameters.output}.locations.bed"
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"

        intersect.biomart_locate(vep_out, self.passParameters)
        self.assertTrue(filecmp.cmp(biomart_out, biomart_pass, shallow=False))

    def test_simulate(self):
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"
        sim_out = f"{self.passParameters.output}.simulation.tmp"
        sim_pass = f"{working_dir}/input/sim_pass.txt"

        binL, binS, bp_ignore, s_ignore = intersect.sim_prep(biomart_pass)
        intersect.simulate(binL, binS, self.passParameters, True)
        self.assertTrue(filecmp.cmp(sim_out, sim_pass, shallow=False))

    def test_simulate_fail(self):
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"
        sim_out = f"{self.passParameters.output}.simulation.tmp"
        sim_pass = f"{working_dir}/input/sim_pass.txt"

        binL, binS, bp_ignore, s_ignore = intersect.sim_prep(biomart_pass)
        intersect.simulate(binL, binS, self.passParameters, False)
        self.assertFalse(filecmp.cmp(sim_out, sim_pass, shallow=False))

    def test_intersect(self):
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"
        l_ignore, s_ignore, bp_total, s_ignore = intersect.sim_prep(biomart_pass)
        int_result = intersect.intersect(self.passParameters, biomart_pass, bp_total)
        self.assertEqual(int_result,(4/2167484))

if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'])#, exit=False)
