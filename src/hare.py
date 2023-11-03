#!/usr/in/env python
# hare.py
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

###############################################################################

# HARE:

# Please cite our article:
# Kun E, Javan EM, Smith OS, et al., The genetic architecture and evolution of the human skeletal form.
# Science 381,eadf8009(2023). DOI:10.1126/science.adf8009

# Dependencies: python => 3.0, pandas, scipy
# See documentation for environment.yml file and links to installation documentation for dependencies

# Inputs: [RESULTS].intersections: Filename for HARE output (can also be comma-separated list) with calculated intersections/bp between HARs for randomly generated controls and provided element set

# Outputs: [output_stem].stats: Tab-separated file with the test parameters and resulting p-value
#          [output_stem].png: Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set

# Example command: hare sigtest --input [INPUT] --output [output_stem]

###############################################################################

import sys
import argparse
import os
from pathlib import Path
from datetime import datetime
now = datetime.now().strftime("%Y%m%d%H%M%S")

def run_intersect(intersect_args):
    import intersect
    intersect.main(**intersect_args)
    return

def run_sigtest(sigtest_args):
    import sigtest
    sigtest.main(**sigtest_args)
    return

def main():
    ########################## RUNNER ARGUMENTS ##########################
    parser = argparse.ArgumentParser(description='Provided GWAS summary statistics, investigate whether there is elevated overlap between genome-wide significant genes and genetic elements of interest.')
    subparsers = parser.add_subparsers(help="HARE functions", dest='command')
    intersect = subparsers.add_parser('intersect', help="Annotate regions, simulate background distribution, and compute intersections/bp between annotated regions and provided genetic elements of interest.")
    sigtest = subparsers.add_parser('sigtest', help="Perform significance testing on results output from intersect function.")

    ########################## INTERSECT ARGUMENTS ##########################
    intersect.add_argument('--gwas', '-g', type=str, help="Complete filepath for GWAS summary statistics.", required=True)
    intersect.add_argument('--eoi', '-e', type=str, help="List of elements of interest (EOIs) in BED format.", required=True, default=None)
    intersect.add_argument('--ref', '-r', type=str, help="Genome reference in BED format. Used to generate random regions for background distribution.", required=True)
    intersect.add_argument('--ref_build', '-b', type=str, help="GRCh reference build. Options are either 37 (hg19) or 38 (hg38). Default is 37.", required=False, default="37")
    intersect.add_argument('--out', '-o', type=str, help="Output filepath stem. Default is the filestem of the input file.", required=False, default=None)
    intersect.add_argument('--pval', '-p', type=float, help="P-value threshold for inclusion of variants in the annotation. Default is 1e-6.", required=False, default=1e-6)
    intersect.add_argument('--maf', '-m', type=float, help="Minimum minor allele frequency (MAF). Default 0.01.", required=False, default=0.01)
    intersect.add_argument('--gwas_p', type=str, help="Column name for p-values in GWAS summary statistics. Default is \"P\".", required=False, default="P")
    intersect.add_argument('--use_z', action="store_true", help="Use Z score (in \"Z\" column) to calculate p-values for GWAS summary statistics. By default is OFF.", required=False)
    intersect.add_argument('--gwas_maf', type=str, help="Column name for minor allele frequency in GWAS summary statistics. Default is \"MAF\".", required=False, default="MAF")
    intersect.add_argument('--gwas_ref', type=str, help="Column name for reference allele in GWAS summary statistics. Default is \"REF\".", required=False, default="REF")
    intersect.add_argument('--gwas_alt', type=str, help="Column name for alternate allele in GWAS summary statistics. Default is \"ALT\".", required=False, default="ALT")
    intersect.add_argument('--snp_map', type=str, help="If SNPs provided as IDs instead of genomic locations (CHR, POS), provide a BED file which maps IDs to locations.", required=False, default=None)
    intersect.add_argument('--source_neale', action="store_true", help="Use Neale Lab GWAS Summary Statistics format (\'GWAS round 2\' version). Will take priority over --gwas_p and --gwas_maf options. Not used by default.")
    intersect.add_argument('--source_bolt', action="store_true", help="Use BOLT-LMM GWAS Summary Statistics format. Will take priority over --gwas_p and --gwas_maf options. Not used by default.")
    intersect.add_argument('--cache_dir', type=str, help="Specify the VEP cache directory to use. Default is \"$HOME/.vep/\"", required=False, default="$HOME/.vep/")
    intersect.add_argument('--cache_version', type=str, help="Specify the VEP cache version to use. Default is 105. It is recommended by Ensembl that you use the same cache and VEP version.", required=False, default="105")
    intersect.add_argument('--dist', '-d', type=int, help="(bp) Distance to transcript for which VEP assigns upstream and downstream consequences. Default is 5000 bp.", required=False, default=5000)
    intersect.add_argument('--biotypes', type=str, help="Allowed biotypes. Options are \"protein_coding\", \"protein_all\", and \"all\". Default is \"protein_all.\" \"all\" will keep all elements regardless of biotype. See VEP's biotype documentation for type details.", required=False, default="protein_all")
    intersect.add_argument('--anno_only', action="store_true", help="Perform only annotation and location finding of genome-wide significant positions (do not simulate and intersect).", required=False)
    # parser.add_argument('--match', type=str, help="Match on recombination rate in addition to length using the provided filepath.", required=False, default=None)
    # parser.add_argument('--non_human', action="store_true", help="Use non-human species. The pipeline will allow this but it is not explicitly supported at this time.")
    intersect.add_argument('--draws', '-n', type=int, help="Number of simulations (draws) to use for background distribution. Default is n=1000.", required=False, default=1000)
    intersect.add_argument('--keep_tmp', action="store_true", help="Keep all intermediate files.")

    ########################## SIGTEST ARGUMENTS ##########################
    sigtest.add_argument('--input','-i', type=str, help="Filepath to results file from HARE with simulation and element set intersection/bp values. Can also be a comma-separated list.", required=True)
    sigtest.add_argument('--out', '-o', type=str, help="Output filepath. Default is \'hare\' if none provided.", required=False, default=None)
    sigtest.add_argument('--distribution', type=str, help="Distribution to use for significance testing. Options are (case insensitive) NORMAL, WEIBULL, BETA, GAMMA. Default is WEIBULL.", required=False, default="weibull")
    sigtest.add_argument('--skip_plot', action="store_true", help="Skip plotting of results. Default is OFF (will plot by default).")

    ########################## RUNNER ##########################
    args = parser.parse_args()
    command = args.command

    if command == "intersect":
        intargs = vars(intersect.parse_known_args()[0])
        run_intersect(intargs)

    elif command == "sigtest":
        sigargs = vars(sigtest.parse_known_args()[0])
        run_sigtest(sigargs)

if __name__ == "__main__":
    main()
