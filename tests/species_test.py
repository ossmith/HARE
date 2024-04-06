#!/usr/in/env python
# species_test.py

# Example command: python species_test.py --cache_dir [VEP_CACHE_PATH] --cache_ver [VEP_VERSION]

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

parser = argparse.ArgumentParser(description='Automated testing of non-human species in HARE pipeline.')
parser.add_argument('--cache_dir', '-c', type=str, help="Specify the VEP cache directory to use. Default is \"$HOME/.vep/\"", required=False, default="$HOME/.vep/")
args = parser.parse_args()

cache_pri = args.cache_dir
cache_dir = os.path.normpath(os.path.expanduser(cache_pri))
# Check directory for VEP cache
if os.path.exists(cache_dir) == False:
    raise FileNotFoundError(f'{cache_dir} does not exist or could not be opened.')
# cache_ver = args.cache_ver
working_dir = os.path.dirname(os.path.realpath(argv[0]))

class TestNonHuman(unittest.TestCase):

    def setUp(self):
        super(TestNonHuman, self).setUp()
        # self.ref = f"{working_dir}/input/ref_pass.bed"
        self.tools = ["vep", "bedtools"]

        # Non-human vertebrate (Bos taurus)
        # Reference cache ARS-UCD1.3 v111: https://ftp.ensembl.org/pub/release-111/variation/vep/bos_taurus_vep_111_ARS-UCD1.3.tar.gz
        self.btaurusRef = f"{working_dir}/input/ref_cow.fai"
        self.btaurusSettings = hareclasses.SettingsContainer(False, False, False, False, True, "bos_taurus", True, False, False, False)
        self.btaurusParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_cow.txt", 1e-6, "P", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test_btaurus", cache_dir, "111", "protein_all", 1000, f"{working_dir}/input/eoi_cow.txt", self.btaurusRef, "37", 3, self.btaurusSettings)

        # Plant (Arabidopsis thaliana)
        # Reference cache TAIR10 v58: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/variation/vep/arabidopsis_thaliana_vep_58_TAIR10.tar.gz
        self.athalianaRef = f"{working_dir}/input/ref_arabidopsis.fai"
        self.athalianaSettings = hareclasses.SettingsContainer(False, False, False, False, True, "arabidopsis_thaliana", False, True, False, False)
        self.athalianaParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_arabidopsis.txt", 1e-6, "p.val", 0.1, "maf", "snp.ref", "snp.alt", None, f"{working_dir}/test_athaliana", cache_dir, "58", "protein_all", 1000, f"{working_dir}/input/eoi_arabidopsis.txt", self.athalianaRef, "37", 3, self.athalianaSettings)

        # Metazoa (Drosophila melanogaster)
        # Reference cache BDGP6.46 v111: https://ftp.ensembl.org/pub/release-111/variation/vep/drosophila_melanogaster_vep_111_BDGP6.46.tar.gz
        self.dmelanogasterRef = f"{working_dir}/input/ref_drosophila.fai"
        self.dmelanogasterSettings = hareclasses.SettingsContainer(False, False, False, False, True, "drosophila_melanogaster", False, False, True, False)
        self.dmelanogasterParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_drosophila.txt", 1e-6, "P-value", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test_dmelanogaster", cache_dir, "111", "protein_all", 1000, f"{working_dir}/input/eoi_drosophila.txt", self.dmelanogasterRef, "37", 3, self.dmelanogasterSettings)

        # Metazoa (Honeybee)
        # Reference cache Amel_HAv3.1 v58: https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/variation/indexed_vep_cache/apis_mellifera_vep_58_Amel_HAv3.1.tar.gz
        self.amelliferaRef = f"{working_dir}/input/ref_bee.fai"
        self.amelliferaSettings = hareclasses.SettingsContainer(False, False, False, False, True, "apis_mellifera", False, False, True, False)
        self.amelliferaParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_bee.txt", 1e-6, "p-value-IND", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test_amellifera", cache_dir, "58", "protein_all", 1000, f"{working_dir}/input/eoi_bee.txt", self.amelliferaRef, "37", 3, self.amelliferaSettings)

        # Fungi (Saccharomyces cerevisiae)
        # Reference cache R64-1-1 v58: https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-58/fungi/variation/vep/saccharomyces_cerevisiae_vep_58_R64-1-1.tar.gz
        self.scerevisiaeRef = f"{working_dir}/input/ref_yeast.fai"
        self.scerevisiaeSettings = hareclasses.SettingsContainer(False, False, False, False, True, "saccharomyces_cerevisiae", False, False, False, True)
        self.scerevisiaeParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_yeast.txt", 1e-6, "P", 0.1, None, "REF", "ALT", None, f"{working_dir}/test_scerevisiae", cache_dir, "58", "protein_all", 1000, f"{working_dir}/input/eoi_yeast.txt", self.scerevisiaeRef, "37", 3, self.scerevisiaeSettings)

        # Bacteria (Staphylococcus aureus)
        # Bacteria and protists not supported by VEP as of version 1.1
        # self.saureusRef = f"{working_dir}/input/ref_staph.bed"
        # self.saureusSettings = hareclasses.SettingsContainer(False, False, False, False, True, "staphylococcus_aureus", False, False, False, False, True)
        # self.saureusParameters = hareclasses.ArgumentContainer(f"{working_dir}/input/gwas_staph.txt", 1e-6, "P", 0.1, "MAF", "REF", "ALT", None, f"{working_dir}/test_staph", cache_dir, 58, "protein_all", 1000, f"{working_dir}/input/eoi_staph.txt", self.saureusRef, "37", 3, self.saureusSettings)

    def tearDown(self):
        # Clean up directory and clear cache
        for f in glob.glob(f"{working_dir}/test_*"):
            os.remove(f)
        filecmp.clear_cache()

    def test_cow(self):
        cow_intpbp = 2/1209043
        cow_filelist = [f"{self.btaurusParameters.output}.snps", f"{self.btaurusParameters.output}.features",
        f"{self.btaurusParameters.output}.locations.bed", f"{self.btaurusParameters.output}.simulation.tmp"]

        # Read and filter GWAS
        cow_snp = intersect.gwas_import(self.btaurusParameters, self.btaurusSettings)
        self.assertTrue(filecmp.cmp(cow_filelist[0], f"{working_dir}/input/cow_pass.snps", shallow=False))

        # VEP annotation
        cow_vep = intersect.vep_annotate(cow_filelist[0], self.btaurusParameters, self.btaurusSettings)
        self.assertTrue(filecmp.cmp(cow_filelist[1], f"{working_dir}/input/cow_pass.features", shallow=False))

        # Find feature locations
        cow_biomart = intersect.biomart_locate(cow_filelist[1], self.btaurusParameters, self.btaurusSettings)
        self.assertTrue(filecmp.cmp(cow_filelist[2], f"{working_dir}/input/cow_pass.locations.bed", shallow=False))

        # Run simulation
        # binL, binS, cow_bp, s_ignore = intersect.sim_prep(f"{working_dir}/input/cow_pass.locations.bed")
        binL, binS, cow_bp, s_ignore = intersect.sim_prep(cow_filelist[2])
        intersect.simulate(binL, binS, self.btaurusParameters, True)
        self.assertTrue(filecmp.cmp(cow_filelist[3], f"{working_dir}/input/cow_pass.sims", shallow=False))

        # Run intersect
        cow_result = intersect.intersect(self.btaurusParameters, cow_filelist[2], cow_bp)
        self.assertEqual(cow_result, cow_intpbp)

    def test_plant(self):
        ara_intpbp = 1/7000
        ara_filelist = [f"{self.athalianaParameters.output}.snps", f"{self.athalianaParameters.output}.features",
        f"{self.athalianaParameters.output}.locations.bed", f"{self.athalianaParameters.output}.simulation.tmp"]

        # Read and filter GWAS
        ara_snp = intersect.gwas_import(self.athalianaParameters, self.athalianaSettings)
        self.assertTrue(filecmp.cmp(ara_filelist[0], f"{working_dir}/input/ara_pass.snps", shallow=False))

        # VEP annotation
        ara_vep = intersect.vep_annotate(ara_filelist[0], self.athalianaParameters, self.athalianaSettings)
        self.assertTrue(filecmp.cmp(ara_filelist[1], f"{working_dir}/input/ara_pass.features", shallow=False))

        # Find feature locations
        ara_biomart = intersect.biomart_locate(ara_filelist[1], self.athalianaParameters, self.athalianaSettings)
        self.assertTrue(filecmp.cmp(ara_filelist[2], f"{working_dir}/input/ara_pass.locations.bed", shallow=False))

        # Run simulation
        binL, binS, ara_bp, s_ignore = intersect.sim_prep(ara_filelist[2])
        intersect.simulate(binL, binS, self.athalianaParameters, True)
        self.assertTrue(filecmp.cmp(ara_filelist[3], f"{working_dir}/input/ara_pass.sims", shallow=False))

        # Run intersect
        ara_result = intersect.intersect(self.athalianaParameters, ara_filelist[2], ara_bp)
        self.assertEqual(ara_result, ara_intpbp)

    def test_drosophila(self):
        dro_intpbp = 1/39154
        dro_filelist = [f"{self.dmelanogasterParameters.output}.snps", f"{self.dmelanogasterParameters.output}.features",
        f"{self.dmelanogasterParameters.output}.locations.bed", f"{self.dmelanogasterParameters.output}.simulation.tmp"]

        # Read and filter GWAS
        dro_snp = intersect.gwas_import(self.dmelanogasterParameters, self.dmelanogasterSettings)
        self.assertTrue(filecmp.cmp(dro_filelist[0], f"{working_dir}/input/dro_pass.snps", shallow=False))

        # VEP annotation
        dro_vep = intersect.vep_annotate(dro_filelist[0], self.dmelanogasterParameters, self.dmelanogasterSettings)
        self.assertTrue(filecmp.cmp(dro_filelist[1], f"{working_dir}/input/dro_pass.features", shallow=False))

        # Find feature locations
        dro_biomart = intersect.biomart_locate(dro_filelist[1], self.dmelanogasterParameters, self.dmelanogasterSettings)
        self.assertTrue(filecmp.cmp(dro_filelist[2], f"{working_dir}/input/dro_pass.locations.bed", shallow=False))

        # Run simulation
        # binL, binS, dro_bp, s_ignore = intersect.sim_prep(f"{working_dir}/input/dro_pass.locations.bed")
        binL, binS, dro_bp, s_ignore = intersect.sim_prep(dro_filelist[2])
        intersect.simulate(binL, binS, self.dmelanogasterParameters, True)
        self.assertTrue(filecmp.cmp(dro_filelist[3], f"{working_dir}/input/dro_pass.sims", shallow=False))

        # Run intersect
        dro_result = intersect.intersect(self.dmelanogasterParameters, dro_filelist[2], dro_bp)
        self.assertEqual(dro_result, dro_intpbp)

    def test_metazoa(self):
        bee_intpbp = 2/1148282
        bee_filelist = [f"{self.amelliferaParameters.output}.snps", f"{self.amelliferaParameters.output}.features",
        f"{self.amelliferaParameters.output}.locations.bed", f"{self.amelliferaParameters.output}.simulation.tmp"]

        # Read and filter GWAS
        bee_snp = intersect.gwas_import(self.amelliferaParameters, self.amelliferaSettings)
        self.assertTrue(filecmp.cmp(bee_filelist[0], f"{working_dir}/input/bee_pass.snps", shallow=False))

        # VEP annotation
        bee_vep = intersect.vep_annotate(bee_filelist[0], self.amelliferaParameters, self.amelliferaSettings)
        self.assertTrue(filecmp.cmp(bee_filelist[1], f"{working_dir}/input/bee_pass.features", shallow=False))

        # Find feature locations
        bee_biomart = intersect.biomart_locate(bee_filelist[1], self.amelliferaParameters, self.amelliferaSettings)
        self.assertTrue(filecmp.cmp(bee_filelist[2], f"{working_dir}/input/bee_pass.locations.bed", shallow=False))

        # Run simulation
        binL, binS, bee_bp, s_ignore = intersect.sim_prep(f"{working_dir}/input/bee_pass.locations.bed")
        binL, binS, bee_bp, s_ignore = intersect.sim_prep(bee_filelist[2])
        intersect.simulate(binL, binS, self.amelliferaParameters, True)
        self.assertTrue(filecmp.cmp(bee_filelist[3], f"{working_dir}/input/bee_pass.sims", shallow=False))

        # Run intersect
        bee_result = intersect.intersect(self.amelliferaParameters, bee_filelist[2], bee_bp)
        self.assertEqual(bee_result, bee_intpbp)

    def test_fungi(self):
        yeast_intpbp = 2/12955
        yeast_filelist = [f"{self.scerevisiaeParameters.output}.snps", f"{self.scerevisiaeParameters.output}.features",
        f"{self.scerevisiaeParameters.output}.locations.bed", f"{self.scerevisiaeParameters.output}.simulation.tmp"]

        # Read and filter GWAS
        yeast_snp = intersect.gwas_import(self.scerevisiaeParameters, self.scerevisiaeSettings)
        self.assertTrue(filecmp.cmp(yeast_filelist[0], f"{working_dir}/input/yeast_pass.snps", shallow=False))

        # VEP annotation
        yeast_vep = intersect.vep_annotate(yeast_filelist[0], self.scerevisiaeParameters, self.scerevisiaeSettings)
        self.assertTrue(filecmp.cmp(yeast_filelist[1], f"{working_dir}/input/yeast_pass.features", shallow=False))

        # Find feature locations
        yeast_biomart = intersect.biomart_locate(yeast_filelist[1], self.scerevisiaeParameters, self.scerevisiaeSettings)
        self.assertTrue(filecmp.cmp(yeast_filelist[2], f"{working_dir}/input/yeast_pass.locations.bed", shallow=False))

        # Run simulation
        binL, binS, yeast_bp, s_ignore = intersect.sim_prep(f"{working_dir}/input/yeast_pass.locations.bed")
        binL, binS, yeast_bp, s_ignore = intersect.sim_prep(yeast_filelist[2])
        intersect.simulate(binL, binS, self.scerevisiaeParameters, True)
        self.assertTrue(filecmp.cmp(yeast_filelist[3], f"{working_dir}/input/yeast_pass.sims", shallow=False))

        # Run intersect
        yeast_result = intersect.intersect(self.scerevisiaeParameters, yeast_filelist[2], yeast_bp)
        self.assertEqual(yeast_result, yeast_intpbp)

    # Bacteria currently not supported by VEP
    # def test_bacteria(self):
    #     staph_filelist = [f"{self.saureusParameters.output}.snps", f"{self.saureusParameters.output}.features",
    #     f"{self.saureusParameters.output}.locations.bed", f"{self.saureusParameters.output}.simulation.tmp"]

if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'])
