#!/usr/in/env python
# HARE.py
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

################################################################################

# This script identifies the genes and genomic features associated with autosomal
# genome-wide significant SNPs based on GWAS summary statistics for a given phenotype.
# It then tests for elevated intersections between the genomic features and
# human accelerated regions (HARs) based on a simulated background distribution.
# This analysis is based upon a modification of the method presented in
# Richard, et al., 2020 (https://doi.org/10.1016/j.cell.2020.02.057).

# Dependencies: python => 3.0, VEP command line tool, human reference genome cache, perl Set::IntervalTree, htslib, bedtools v2.30.0

# Inputs: Files containing summary statistics and references for intersection testing.
#         - [GWAS]: GWAS summary statistics file which contains information on the chromosome, position, and p-value for each SNP
#         - [HAR].bed: BED file with human accelerated regions (HARs)
#         - [REF].bed: BED file with gene annotation of human genome reference (e.g. UCSC's hg19 .bed)

# Outputs: File set in HARE_results directory.
#       - [OUT_STEM].snps: variants which pass QC conditions (biallelic, MAF threshold, p-value threshold)
#       - [OUT_STEM].annotation: raw VEP annotation results
#       - [OUT_STEM].annotation_summary.html: VEP run log produced by command line run
#       - [OUT_STEM].biomart: Output from Biomart query containing Ensembl gene ID, gene start and end position, chromosome, gene name (symbol), and strand.
#       - [OUT_STEM].locations.bed: BED file containing all the locations of the genes selected for inclusion in the element set. Input for har_test.py script to compute HAR intersections.
#       - [OUT_STEM].intersections: calculated intersections/bp between HARs for randomly generated controls and provided element set

# Example command: python HARE.py --gwas [GWAS] --HAR [HAR] --ref [REF] --out [OUT_STEM] [OPTIONS]...

################################################################################

#### Import dependencies and setup ####
import argparse
import pandas as pd
import os
from os import path
import subprocess
from datetime import datetime
import numpy as np
import ntpath
import random
from io import StringIO

script_version = "1.0.0"

parser = argparse.ArgumentParser(description='Provided GWAS summary statistics, investigate whether there is elevated overlap between genome-wide significant genes and human accelerated regions (HARs).')
parser.add_argument('--gwas', '-g', type=str, help="Complete filepath for GWAS summary statistics. FILES MUST BE TAB SEPARATED.", required=True)
parser.add_argument('--out', '-o', type=str, help="Output filepath stem. Default is the filestem of the GWAS summary statistics file.", required=False, default=None)
parser.add_argument('--ref', '-r', type=str, help="Genome reference in BED format. Used to generate random regions for background distribution.", required=True)
parser.add_argument('--ref_build', '-b', type=str, help="GRCh reference build. Options are either 37 (hg19) or 38 (hg38).", required=False, default="37")
parser.add_argument('--pval',"-p", type=float, help="P-value threshold for inclusion of variants in the annotation. Default is 1e-6.", required=False, default=1e-6)
parser.add_argument('--maf', type=float, help="Minimum minor allele frequency (MAF). Default 0.01.", required=False, default=0.01)
parser.add_argument('--gwas_p', type=str, help="Column name for p-values in GWAS summary statistics. Default is \"P\".", required=False, default="P")
parser.add_argument('--gwas_maf', type=str, help="Column name for minor allele frequency (MAF) in GWAS summary statistics. Default is \"MAF\".", required=False, default="MAF")
parser.add_argument('--source_neale', type=str, help="Use Neale Lab GWAS Summary Statistics format. Will take priority over --gwas_p and --gwas_maf options. Accepts T, True, F, and False (case insensitive). Default False.", required=False, default="F")
parser.add_argument('--source_bolt', type=str, help="Use BOLT-LMM GWAS Summary Statistics format. Will take priority over --gwas_p and --gwas_maf options. Accepts T, True, F, and False (case insensitive). Default False.", required=False, default="F")
parser.add_argument('--cache_dir', type=str, help="Specify the VEP cache directory to use. Default is \"$HOME/.vep/\"", required=False, default="$HOME/.vep/")
parser.add_argument('--cache_version', type=str, help="Specify the VEP cache version to use. Default is 105 (version used during development).", required=False, default="105")
parser.add_argument('--dist', "-d", type=int, help="(bp) Distance to transcript for which VEP assigns upstream and downstream consequences. Default is 5,000 bp.", required=False, default=5000)
parser.add_argument('--biotypes', type=str, help="Allowed biotypes. Options are \"protein_coding\", \"protein_all\", and \"all_features\". Default is \"protein_all.\" See VEP's biotype documentation for details.", required=False, default="protein_all")
parser.add_argument('--HAR',"-H", type=str, help="List of human accelerated regions (HARs) in BED format. Default is \'harsRichard2020.bed\' (included with GitHub repo.)", required=False, default="harsRichard2020.bed")
parser.add_argument('--draws', '-n', type=int, help="Number of simulations (draws) to use for background distribution. Default is n=1,000.", required=False, default=1000)
args = parser.parse_args()
NaN = np.nan

# Arguments for VEP and BioMart annotation
gwas_filepath = args.gwas
if os.path.exists(gwas_filepath) == False:
    print(f'\nFileNotFoundError: {gwas_filepath} does not exist or could not be opened.')
    exit()
pval = args.pval
pval_col = args.gwas_p
maf = args.maf
maf_col = args.gwas_maf
output_stem = args.out
if output_stem == None:
    output_stem = os.path.splitext(ntpath.basename(gwas_filepath))[0]
cache_dir = args.cache_dir
if os.path.exists(cache_dir) == False:
    print(f'\nFileNotFoundError: {cache_dir} does not exist or could not be opened.')
    exit()
cache_v = args.cache_version
biotypes = args.biotypes
dist = args.dist
source_neale = args.source_neale[0].upper() # Uses capitalized first letter to determine true or false
source_bolt = args.source_bolt[0].upper()
if (source_neale == "T" and source_bolt == "T"):
    print("Note: source_neale and source_bolt cannot both be set to True. Defaulting to Neale Lab format.")
if source_neale == "T":
    pval_col = "pval"
    maf_col = "minor_AF"
elif source_bolt == "T":
    pval_col = "P_BOLT_LMM_INF"
    maf_col = "MAF"

# Arguments for HAR testing
har = args.HAR
if os.path.exists(har) == False:
    print(f'\nFileNotFoundError: {har} does not exist or could not be opened.')
    exit()
ref = args.ref
if os.path.exists(ref) == False:
    print(f'\nFileNotFoundError: {ref} does not exist or could not be opened.')
    exit()
build = args.ref_build
if build not in ["37","38"]:
    print("[input] Genome build must be either 37 or 38! Using default of GRCh37.")
draws = args.draws

from pathlib import Path
RESULTS_DIR = f"HARE_results"
Path(RESULTS_DIR).mkdir(parents=True, exist_ok=True)
output = f"{RESULTS_DIR}/{output_stem}"

#### Import and QC the GWAS summary statistics to identify genome-wide significant variants ####
def gwas_import(gwas):
    print('Loading GWAS summary statistics...', end="", flush=True)
    if source_neale == "T":
        raw_df = pd.read_csv(gwas, sep="\t", header=0, index_col=False)
        raw_df[["CHROM","POS","REF","ALT"]] = raw_df["variant"].str.split(":",expand=True)
    elif source_bolt == "T":
        # headers = ["ID","CHROM","POS","GENPOS","ALT","REF","A1FREQ","F_MISS","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM"]
        raw_df = pd.read_csv(gwas, sep="\t", header=0, usecols=["CHR","BP","SNP","ALLELE1","ALLELE0","P_BOLT_LMM_INF"], index_col=False)
        raw_df = raw_df.rename(columns={"CHR":"CHROM", "BP":"POS", "SNP":"ID", "ALLELE1":"ALT", "ALLELE0":"REF"})
    else:
        raw_df = pd.read_csv(gwas, sep="\t", header=0, index_col=False)
    print("OK")

    print("Filtering SNPs for analysis...",end="", flush=True)
    raw_df.columns = raw_df.columns.str.replace('[#,@,&]', '',regex=True)

    gwas_df = raw_df[(raw_df["CHROM"]!="X") & (raw_df["CHROM"]!="Y") & (raw_df["CHROM"]!="M") & (raw_df["CHROM"]!=23) & (raw_df["CHROM"]!="23")] # Autosomes only
    pd.options.mode.chained_assignment = None # Because this chained assignment is ok, we will suppress this warning
    gwas_df[pval_col] = gwas_df[pval_col].replace(to_replace=list([NaN,"-"]),value=1)
    gwas_df = gwas_df[gwas_df[pval_col]<pval] # p-val threshold

    if maf_col in gwas_df.columns: # MAF threshold
        gwas_df = gwas_df[(gwas_df[maf_col]>maf)]
    if ("REF" in gwas_df.columns) & ("ALT" in gwas_df.columns): # Biallelic only
        gwas_df = gwas_df[((gwas_df["REF"].isin(["A","T","C","G"])) & (gwas_df["ALT"].isin(["A","T","C","G"])))]
    last_col = list(["ID","QUAL","FILTER","INFO","REF","ALT"])
    for c in last_col:
        if c not in gwas_df.columns:
            gwas_df[c] = "."
        else:
            gwas_df[c] = gwas_df[c].fillna(".")

    gwas_df = gwas_df.sort_values(axis=0, by=["CHROM","POS"], ascending=True)
    snps_out = f"{output}.snps"
    gwas_df[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]].to_csv(snps_out, sep="\t", index=False, header=False)

    print("OK")
    return snps_out


#### Use Ensembl GRCh37 Variant Effect Predictor to identify features which are within {dist} (default 5,000 bp) of genome-wide significant SNPs ####
def vep_annotate(snps_out):
    print("Annotating SNPs...", end="", flush=True)
    vep_out = f"{output}.annotation"
    ann_names = ["VAR_PROVIDED", "LOCATION", "ALLELE", "GENE", "FEATURE", "FEATURE_TYPE", "CONSEQUENCE", "cDNA_POS", "CDS_POS", "PROT_POS", "AMINO_ACIDS", "CODONS", "EXISTING_VARIATION", "EXTRA"]
    if biotypes == "regulatory":
        cmd_annotate = ['vep', '--offline', '--cache', '--dir_cache', cache_dir, '-i', snps_out, '-o', vep_out, '--cache_version', cache_v, '--distance', str(dist), '--biotype','--force_overwrite','--regulatory']
    else:
        cmd_annotate = ['vep', '--offline', '--cache', '--dir_cache', cache_dir, '-i', snps_out, '-o', vep_out, '--nearest', 'gene', '--cache_version', cache_v, '--distance', str(dist), '--biotype','--force_overwrite']
    run_annotate = subprocess.run(cmd_annotate, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Remove duplicate genes and limit to only genes in the desired biotype
    try:
        annotation_df = pd.read_csv(vep_out, sep="\t", skipinitialspace=True, index_col=False, names=ann_names, usecols=["LOCATION","ALLELE","GENE","FEATURE","FEATURE_TYPE","EXTRA"], comment='#')
    # print(annotation_df)
    except FileNotFoundError:
        cmd_msg = " ".join(cmd_annotate)
        print(f'\nFileNotFoundError: VEP annotation output file not found. This often occurs when there are issues with the VEP cache or version compatibility. To troubleshoot, try running the command: {cmd_msg}. Exiting.')
        exit()

    annotation_df[["X","BIOTYPE"]] = annotation_df["EXTRA"].str.split("BIOTYPE=", expand=True)
    if biotypes == "protein_coding" or biotypes == "protein_all":
        annotation_df[["BIOTYPE","NEAREST"]] = annotation_df["BIOTYPE"].str.split(";NEAREST=", expand=True)
    if biotypes == "protein_coding":
        allowed_biotypes=list(["protein_coding"])
        annotation_df = annotation_df[annotation_df["BIOTYPE"].isin(allowed_biotypes)]
        dup_drop = "NEAREST"
    elif biotypes == "protein_all":
        allowed_biotypes=list(["protein_coding","IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_M_gene","IG_V_gene","IG_Z_gene","nonsense_mediated_decay","nontranslating_CDS","non_stop_decay","polymorphic_pseudogene","TR_C_gene","TR_D_gene","TR_J_gene"])
        annotation_df = annotation_df[annotation_df["BIOTYPE"].isin(allowed_biotypes)]
        dup_drop = "NEAREST"
    elif biotypes == "regulatory":
        annotation_df = annotation_df[annotation_df["FEATURE_TYPE"]=="RegulatoryFeature"]
        dup_drop = "FEATURE"
        print(annotation_df)
    annotation_df = annotation_df.drop_duplicates(subset=dup_drop).drop(columns=["X","EXTRA"])

    annotation_out = f"{output}.features"
    annotation_df.to_csv(annotation_out, sep="\t", index=False, header=True)
    print("OK")
    return annotation_out

#### Retrieve genome location (CHR, START, END) for each element in the annotation list ####
def biomart_locate(annotation_out):
    print("Locating elements with BioMart...",end="", flush=True)

    if biotypes == "regulatory":
        feature_col = "FEATURE"
        biomart_header = ["CHR","START","END","FEATURE_TYPE","FEATURE_ID"]
    else:
        feature_col = "NEAREST"
        biomart_header = ["ENSEMBL_ID","START","END","CHR","GENE_NAME","STRAND"]

    # Grab the location (CHR and POS) of the genes identified through vep_annotate()
    annotation_ids = pd.read_csv(annotation_out, sep="\t", header=0, usecols=[feature_col])
    biomart_out = f"{output}.biomart"
    biomart_df = pd.DataFrame(columns=biomart_header)
    biomart_tmp = f"{biomart_out}.tmp"
    ids_count = len(annotation_ids[feature_col])

    if build == "38":
        ensembl_build = "ensembl.org"
    else:
        ensembl_build = "grch37.ensembl.org"

    # Biomart doesn't like to accept lists longer than 400, so we will break it up if this is the case
    if ids_count<=400:
        genes = list(annotation_ids[feature_col])
        ensembl_ids = ",".join(genes)

        if biotypes == "regulatory":
            cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_regulatory_feature\" interface = \"default\" ><Filter name = \"regulatory_stable_id\" value = \"{ensembl_ids}\"/><Attribute name = \"chromosome_name\" /><Attribute name = \"chromosome_start\" /><Attribute name = \"chromosome_end\" /><Attribute name = \"feature_type_name\" /><Attribute name = \"regulatory_stable_id\" /></Dataset></Query>\'"
        else:
            cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"ensembl_gene_id\" value = \"{ensembl_ids} \"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"strand\" /></Dataset></Query>\'"
        os.system(cmd_biomart)
        tmp_df = pd.read_csv(biomart_tmp, sep="\t", names=biomart_header)
        biomart_df = biomart_df.append(tmp_df)
    else:
        range_list = list(range(400,ids_count,400))
        range_list.insert(0,0)
        if range_list[-1] != ids_count:
            range_list.append(ids_count)
        for i in range(len(range_list)-1):
            start = range_list[i]
            end = range_list[i+1]
            genes = list(annotation_ids[feature_col])[start:end]
            ensembl_ids = ",".join(genes)
            if biotypes == "regulatory":
                cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_regulatory_feature\" interface = \"default\" ><Filter name = \"regulatory_stable_id\" value = \"{ensembl_ids}\"/><Attribute name = \"chromosome_name\" /><Attribute name = \"chromosome_start\" /><Attribute name = \"chromosome_end\" /><Attribute name = \"feature_type_name\" /><Attribute name = \"regulatory_stable_id\" /></Dataset></Query>\'"
            else:
                cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"ensembl_gene_id\" value = \"{ensembl_ids} \"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"strand\" /></Dataset></Query>\'"
            os.system(cmd_biomart)
            tmp_df = pd.read_csv(biomart_tmp, sep="\t", names=biomart_header)
            biomart_df = biomart_df.append(tmp_df)
    biomart_df.to_csv(biomart_out, sep="\t", header=True, index=False)

    print("OK")
    print(f"[biomart] Total number of annotations: {ids_count}")

    # Reformat to bed file for import into bedtools and HAR enrichment testing
    locations_out = f"{output}.locations.bed"
    print(f"Writing output to {locations_out}...", end="", flush=True)
    locations_df = pd.DataFrame({})
    if build == "38":
        biomart_df["CHR"] = "chr" + biomart_df["CHR"].astype(str)
    locations_df[1], locations_df[2], locations_df[3] = biomart_df["CHR"], biomart_df["START"], biomart_df["END"]
    locations_df.to_csv(locations_out, sep="\t", header=False, index=False)
    print("OK")
    return locations_out

def read_in(locations_out):
    print("Loading element set...",end="", flush=True)
    df_elemSet = pd.read_csv(locations_out, skipinitialspace=True, sep="\t", header=None, names=["CHR","START","END"], usecols=[0,1,2])
    print("OK")
    sSize = len(df_elemSet.index) # Get the size of the element set
    eDiff = list()
    df_elemSet["LENGTH"] = abs(df_elemSet["END"]-df_elemSet["START"])
    df_elemSet["BIN"] = pd.qcut(df_elemSet["LENGTH"],q=10,duplicates="drop")
    bin_len = df_elemSet.groupby(["BIN"]).mean()["LENGTH"]
    bin_size = df_elemSet.groupby(["BIN"]).count()["LENGTH"] # All columns have the same value which is a count of how many elements are in the bin
    bin_bp = df_elemSet["LENGTH"].sum()
    lens = list()
    sizes = list()
    for i in range(len(bin_len)):
        lens.append(np.floor(bin_len[i]))
        sizes.append(bin_size[i])
    print("[elements] Average element lengths: ", lens)
    print("[elements] Number of elements in each bin: ", sizes)
    print(f"[elements] Element set size: {bin_bp} bp")
    return bin_len, bin_size, bin_bp, sSize

#### Generate background distribution by simulating matched element sets ####
def simulate(bin_len, bin_size, refGenome):
    simulation_outPath = f"{output}.simulation.tmp"

    with open(simulation_outPath, "w") as outfile:
        for i in range(0,len(bin_len)):
          L = str(bin_len[i])
          S = str(bin_size[i])
          cmd_rand = ["bedtools","random","-l",L,"-n",S,"-g",refGenome]
          run_rand = subprocess.run(cmd_rand, stdout=outfile, stderr=subprocess.PIPE)

    outfile.close()
    return simulation_outPath

#### Find intersections between element set and HARs ####
def intersect(fileHAR,fileB,eBP):
    intersect_outPath = f"{output}.intersect.tmp"
    cmd_int = ["bedtools","intersect","-c","-a",str(fileB),"-b",str(fileHAR)]
    with open(intersect_outPath, "w") as outfile:
        run_int=subprocess.Popen(cmd_int,stdin=subprocess.PIPE,stdout=outfile)
        stdout, stderr = run_int.communicate()

    with open(intersect_outPath,"r") as f:
        int_count = 0
        for line in f:
            int_count += int(line.split()[-1])

    int_per_bp = int_count/eBP
    return int_per_bp

def main():
    print("\n ----------------------------------------------------------------------")
    print("|                                                                      |")
    print("|    \U0001F430 Welcome to the HAR Enrichment Analysis Pipeline (HARE) \U0001F430.     |")
    print("|  Please contact Olivia Smith at osmith@utexas.edu to report issues.  |")
    print("|                                                                      |")
    print(" ----------------------------------------------------------------------\n")

    print(f"\U0001F955-----Workflow started at {datetime.now()}.\n")
    print(f"[args] GWAS summary statistics file: \'{gwas_filepath}\'")
    if source_neale == "T":
        print(f"[args] Using Neale Lab summary statistics format.")
    elif source_bolt == "T":
        print(f"[args] Using BOLT-LMM summary statistics format.")
    print(f"[args] P-value significance threshold: {pval}")
    print(f"[args] Minor allele frequency (MAF) threshold: {maf}")
    print(f"[args] VEP cache directory: \'{cache_dir}\', version {cache_v}")
    print(f"[args] Annotation max distance: {dist}")
    print(f"[args] Biotypes included in analysis: {biotypes}")
    print(f"[args] Reference genome annotation file: \'{ref}\'")
    print(f"[args] Reference genome build: GRCh{build}")
    print(f"[args] Human accelerated regions (HAR) BED file: \'{har}\'")
    print(f"[args] Simulation draws to be performed: {draws}")

    # Annotation
    snps_filepath = gwas_import(gwas_filepath)
    annotation_filepath = vep_annotate(snps_filepath)
    locations_filepath = biomart_locate(annotation_filepath)

    # HAR intersection testing
    intersections = list()
    cat = list()
    len_vector, s_vector, total_bp, set_size = read_in(locations_filepath)

    # Generate and test simulation set
    print("Simulating matched element sets and identifying intersections...", end="", flush=True)
    for d in range(draws):
        simulation_results = simulate(len_vector, s_vector, ref)
        intersections.append(intersect(har, simulation_results, total_bp))
        cat.append("simulation")
    print("OK")

    # Test element set of interest
    print("Examining test element set...",end="", flush=True)
    set_intersect = intersect(har, locations_filepath, total_bp)
    intersections.append(set_intersect)
    cat.append("test_set")
    print("OK")

    #### Generate output file ####
    print("Generating output file...",end="", flush=True)
    out_df = pd.DataFrame({"category":cat, "int_per_bp":intersections, "set_size":set_size})
    out_path = f"{output}.intersections"
    out_df.to_csv(out_path,sep="\t",index=False)
    print("OK")

    ### Remove temporary files ###
    print("Cleaning up directory...",end="", flush=True)
    os.remove(f"{output}.biomart.tmp")
    os.remove(simulation_results)
    os.remove(f"{output}.intersect.tmp")
    print("OK")

    print(f"\nWorkflow completed at {datetime.now()}.-----\U0001F407\n")

if __name__ == "__main__":
    main()
