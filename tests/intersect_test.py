#!/usr/in/env python
# intersect_test.py

import unittest
import intersect
import os
import filecmp
import hareclasses
import argparse
from sys import exit
from sys import argv
import glob

parser = argparse.ArgumentParser(description='Automated testing of HARE pipeline.') # Must allow specification of VEP cache
parser.add_argument('--cache_dir', '-c', type=str, help="Specify the VEP cache directory to use. Default is \"$HOME/.vep/\"", required=False, default="$HOME/.vep/")
parser.add_argument('--cache_ver', '-v', type=str, help="Specify the VEP cache version to use. Default is 105. It is recommended by Ensembl that you use the same cache and VEP version.", required=False, default="105")
args = parser.parse_args()

cache_pri = args.cache_dir
cache_dir = os.path.normpath(os.path.expanduser(cache_pri))
cache_ver = args.cache_ver
working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestAnnotation(unittest.TestCase):

    def __init__(self, *arguments):
        super().__init__()

    def runTest(self):
        ref = f"{working_dir}/input/ref_pass.bed" # Going to use the same annotation reference file for all

        # Pass
        passSettings = hareclasses.SettingsContainer(None, None, None, None, None)
        passParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_pass.txt", 1e-6, "P", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", ref, "37", 3, passSettings)

        # Neale
        nealeSettings = hareclasses.SettingsContainer(None, None, True, None, None)
        nealeParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_neale.txt", 1e-6, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", ref, "37", 3, nealeSettings)

        # Bolt
        boltSettings = hareclasses.SettingsContainer(None, None, None, True, None)
        boltParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_bolt.txt", 1e-6, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_pass.txt", ref, "37", 3, boltSettings)

        # Fail
        failParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_fail.txt", 1e-100, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test", cache_dir, cache_ver, "protein_all", 1000, f"{working_dir}/input/eoi_fail.txt", ref, "38", 3, passSettings)

        allParameters = [passParameters, nealeParameters, boltParameters]
        allSettings = [passSettings, nealeSettings, boltSettings]

        snpFiles = [f"{working_dir}/input/snp_pass.snps", f"{working_dir}/input/snp_neale.snps", f"{working_dir}/input/snp_bolt.snps"]

        for fileset in allParameters: # Make sure the gwas and ref files exist
            self.find_files(fileset)

        self.test_header(allSettings, allParameters)
        self.test_header_fail(passSettings, failParameters)
        # self.test_snp2loc(testSettings, passParameters)
        self.test_gwas(allSettings, allParameters, snpFiles)
        self.test_gwas_fail(failParameters, passSettings)
        self.test_vep(passParameters)
        self.test_biomart(passParameters)
        self.test_simulate(passParameters)
        self.test_intersect(passParameters)
        self.clean_dir() # Clean directory after
        pass

    def find_files(self, params):
        for f in [params.gwas, params.cache_dir, params.eoi, params.ref]:
            intersect.findFiles(f)
        return

    # def test_snp2loc(self, testSettings, passParameters):
    #     CHECK SNP MAP
    #     self.assertEqual(intersect.snpToLoc(testSettings, passParameters))
    #
    #     with self.assertRaises(RuntimeError):
    #         intersect.snpToLoc(self) # Mapping failed because mapped_df = 0
    #
    #     return

    def test_header(self, settingsSet, paramSet):
        defaultHeader = ["CHR", "POS", "P"]
        nealeHeader = ["variant", "pval"]
        boltHeader = ["BP", "SNP", "ALLELE1", "ALLELE0", "P_BOLT_LMM"]
        headersList = [defaultHeader, nealeHeader, boltHeader]
        colList = [["CHR","POS"], ["CHR","POS"], ["CHR","BP"]]

        for h in range(len(headersList)):
            chrResult, posResult = intersect.check_header(paramSet[h], settingsSet[h], headersList[h], colList[h][0], colList[h][1])
            self.assertEqual(chrResult, colList[h][0])
            self.assertEqual(posResult, colList[h][1])
        return

    def test_header_fail(self, testSettings, failParameters):

        with self.assertRaises(KeyError):
            intersect.check_header(failParameters, testSettings, ["CHR", "POS", "FOO"], "CHR", "POS")

        return

    def test_gwas(self, settingsClass, paramSet, snpSet):

        for g in range(len(snpSet)):
            snpLoc = intersect.gwas_import(paramSet[g], settingsClass[g])
            self.assertTrue(filecmp.cmp(snpLoc, snpSet[g], shallow=False))

        return

    def test_gwas_fail(self, argumentClass, settingsClass):

        with self.assertRaises(RuntimeError):
            intersect.gwas_import(argumentClass, settingsClass)

        return

    def test_vep(self, argumentClass):
        snps_out = f"{working_dir}/input/snp_pass.snps"
        vep_out = f"{argumentClass.output}.features"
        vep_success = f"{working_dir}/input/vep_pass.features"

        intersect.vep_annotate(snps_out, argumentClass)
        self.assertTrue(filecmp.cmp(vep_out, vep_success, shallow=False))

        return

    def test_biomart(self, argumentClass):
        vep_out = f"{working_dir}/input/vep_pass.features"
        biomart_out = f"{argumentClass.output}.locations.bed"
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"

        intersect.biomart_locate(vep_out, argumentClass)
        self.assertTrue(filecmp.cmp(biomart_out, biomart_pass, shallow=False))

        return

    def test_simulate(self, argumentClass):
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"
        sim_out = f"{argumentClass.output}.simulation.tmp"
        sim_pass = f"{working_dir}/input/sim_pass.txt"

        binL, binS, bp_ignore, s_ignore = intersect.sim_prep(biomart_pass)
        intersect.simulate(binL, binS, argumentClass, True)
        self.assertTrue(filecmp.cmp(sim_out, sim_pass, shallow=False))

        return

    def test_intersect(self, argumentClass):
        biomart_pass = f"{working_dir}/input/biomart_pass.locations.bed"
        l_ignore, s_ignore, bp_total, s_ignore = intersect.sim_prep(biomart_pass)
        int_result = intersect.intersect(argumentClass, biomart_pass, bp_total)
        self.assertEqual(int_result,(4/2167484))
        # pass
        return

    def clean_dir(self):
        for f in glob.glob(f"{working_dir}/test.*"):
            os.remove(f)
        filecmp.clear_cache() # Remove filecmp cache
        return

if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'])#, exit=False)
