#!/usr/in/env python
# intersect.py
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

###############################################################################

# This script identifies the genes and genomic features associated with autosomal
# genome-wide significant SNPs based on GWAS summary statistics for a given phenotype.
# It then tests for elevated intersections between the genomic features and
# elements of interest (EOIs) based on a simulated background distribution.
# This analysis is an expansion of the method presented in
# Richard, et al., 2020 (DOI:10.1016/j.cell.2020.02.057).

# Dependencies: python => 3.0, VEP command line tool, human reference genome cache, perl Set::IntervalTree, htslib, bedtools v2.30.0, wget
# See documentation for environment.yml file and links to installation documentation for dependencies

# Inputs: Files containing summary statistics and references for intersection testing.
#         - [GWAS]: GWAS summary statistics file which contains information on the chromosome, position, and p-value for each SNP
#         - [EOI].bed: BED file with elements of interest (EOIs)
#         - [REF].bed: BED file with gene annotation of human genome reference (e.g. UCSC's hg19 .bed)

# Outputs: File set in HARE_results directory.
#       - [OUT_STEM].snps: variants which pass QC conditions (biallelic, MAF threshold, p-value threshold)
#       - [OUT_STEM].annotation: raw VEP annotation results
#       - [OUT_STEM].annotation_summary.html: VEP run log produced by command line run
#       - [OUT_STEM].biomart: Output from Biomart query containing Ensembl gene ID, gene start and end position, chromosome, gene name (symbol), and strand.
#       - [OUT_STEM].locations.bed: BED file containing all the locations of the genes selected for inclusion in the element set. Input for har_test.py script to compute HAR intersections.
#       - [OUT_STEM].intersections: calculated intersections/bp between element set for randomly generated controls and provided element set

# Example command: hare intersect --gwas [GWAS] --eoi [EOI] --ref [REF] --out [OUT_STEM] [OPTIONS]...

###############################################################################

###############################################################################
################ Setup libraries, arguments, and dependencies #################
###############################################################################
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
from shutil import which
from pathlib import Path
import scipy.stats as stats
import re
import hareclasses
NaN = np.nan

###############################################################################
################# GWAS import and preparation for annotation ##################
###############################################################################
def findFiles(filename):
    if filename == None:
        return
    if os.path.exists(filename) == False:
        raise FileNotFoundError(f'{filename} does not exist or could not be opened.')
    return

def snpToLoc(gwas_df, map):
    '''
    rsIDs are not allowed when using the offline VEP command line interface, so
    must use a map to convert rsID to CHR:POS.
    '''
    snp_df = pd.read_table(map, header=0, comment="##")
    # snp_df = snp_df.rename(columns={0:"CHR", 1:"POS", 2:"END", 3:"SNP"})
    mapped_df = snp_df[["CHR","POS","SNP"]].merge(gwas_df, how="inner", on="SNP").drop_duplicates(subset=["SNP"], ignore_index=True)

    # Check that all SNPs are mapped properly, if not return warning
    if len(mapped_df) == 0:
        print(f"\nRuntimeError: Mapping failed. Confirm that your map BED file contains the SNPs present in the GWAS. \nExiting.\n")
        exit()
    if len(gwas_df) != len(mapped_df):
        lostSNPs = len(gwas_df) - len(mapped_df)
        print(f"\nWARNING: Not all SNP IDs successfully mapped to genomic coordinates. {lostSNPs} were not mapped and will be removed.")
    # Remove 'chr' prefix for GRCh37 build
    if build == "37":
        mapped_df["CHR"] = mapped_df["CHR"].apply(lambda s: s.replace("chr", ""))
    # mapped_df.to_csv("TESTINGMAP.tsv",sep="\t",header=True,index=False)
    return mapped_df

def check_header(argumentClass, settingsClass, headerList, chr_col, pos_col):#, not_human)
    '''
    Confirm that the GWAS summary statistics and EOI files have an acceptable
    format prior to loading them. The formatting and version matching is
    critical to ensure the various tools run/return correct results.
    '''

    dfHead = pd.read_csv(argumentClass.gwas, header=0, sep="\t", index_col=False, nrows=1)

    if (settingsClass.source_neale == False) & (settingsClass.source_bolt == False) & (argumentClass.snp_map == None):
        chr_col, pos_col = "CHR", "POS"
        # Allow CHR or CHROM and BP or POS
        if ("CHR" not in dfHead.columns) & ("CHROM" in dfHead.columns):
            headerList[headerList.index("CHR")] = "CHROM"
            chr_col = "CHROM"
        if ("POS" not in dfHead.columns) & ("BP" in dfHead.columns):
            headerList[headerList.index("POS")] = "BP"
            pos_col = "BP"

    # print(dfHead.columns)
    # print(headerList)
    #
    # if all(item in dfHead.columns for item in headerList) == True:
    #     print("it worked")

    if all(item in dfHead.columns for item in headerList) == False:
        missing = list(set(headerList) - set(dfHead.columns))
        raise KeyError(f"Necessary columns are not present in GWAS file. Columns specifying chromosome, position, reference allele, and alternate allele are all required. \n{missing} columns are missing \nSpecify different column headers using --gwas_<HEADER>. Also note that the file must be TAB (\\t) separated. \nExiting.")
        # print("\nKeyError: Necessary columns are not present in GWAS file. Columns specifying chromosome, position, reference allele, and alternate allele are all required.")
        # print(f"{missing} columns are missing.")
        # print("Specify different column headers using --gwas_{HEADER}. Also note that the file must be TAB (\\t) separated. \nExiting.")
        # exit()

    if settingsClass.anno_only == False:
        # Check that the chromosome naming is consistent with that of the specified GRCh build.
        # Only read a small number of lines because we want to check the file formatting, not read it into memory
        eoidf = pd.read_csv(argumentClass.eoi, delim_whitespace=True, nrows=50)
        try: # Using try/except because will fail if that column isn't a string, we can allow that and set eoi_check to 0
            eoi_check = eoidf.iloc[:,0].str.contains('chr').sum()
        except:
            eoi_check = 0
        #if not_human == True:
            #return
        if (argumentClass.build == 38) & (eoi_check == 0):
            print("\nSyntaxError: Inconsistent chromosome naming for GRCh38 in EOI file. Chromosome names must use \'chr\' prefix. Please check that you are using the correct reference build. Exiting.")
            exit()
        elif (argumentClass.build == 37) & (eoi_check != 0):
            print("\nSyntaxError: Inconsistent chromosome naming for GRCh37 in EOI file. Chromosome names must not use \'chr\' prefix. Please check that you are using the correct reference build. Exiting.")
            exit()

    return chr_col, pos_col

def gwas_import(argumentClass, settingsClass):
    '''
    Read in GWAS summary statistics file, map rsID to CHR:POS if necessary, and
    use given conditions (biallelic, autosome, MAF and p-value thresholds) to
    select which SNPs to annotate.
    '''
    print('Loading GWAS summary statistics...', end="", flush=True)
    if settingsClass.source_neale == True:
        nealeHeader = ["variant", argumentClass.p_col, argumentClass.maf_col]
        check_header(argumentClass, settingsClass, nealeHeader, None, None)
        raw_df = pd.read_csv(argumentClass.gwas, sep="\t", header=0, index_col=False)
        raw_df[["CHR","POS","REF","ALT"]] = raw_df["variant"].str.split(":",expand=True)
    elif settingsClass.source_bolt == True:
        boltHeader = ["BP", "SNP", "ALLELE1", "ALLELE0", argumentClass.p_col] #, maf_col]
        check_header(argumentClass, settingsClass, boltHeader, "CHR", "BP")
        raw_df = pd.read_csv(argumentClass.gwas, sep="\t", header=0, usecols=["CHR","BP","SNP","ALLELE1","ALLELE0",argumentClass.p_col], index_col=False)
        raw_df = raw_df.rename(columns={"BP":"POS", "SNP":"ID", "ALLELE1":"ALT", "ALLELE0":"REF"})
    else:
        defaultHeader = ["CHR", "POS", argumentClass.p_col, argumentClass.alt_col, argumentClass.ref_col]
        snpHeader = ["SNP", argumentClass.p_col, argumentClass.alt_col, argumentClass.ref_col]
        if argumentClass.snp_map == None:
            chr_col, pos_col = check_header(argumentClass, settingsClass, defaultHeader, "CHR", "POS")
        else:
            chr_col, pos_col = check_header(argumentClass, settingsClass, snpHeader, "CHR", "POS")
        raw_df = pd.read_csv(argumentClass.gwas, sep="\t", header=0, index_col=False) # sep="\t"
        if chr_col != "CHR":
            raw_df.rename(columns={chr_col:"CHR"}, inplace=True)
        if pos_col != "POS":
            raw_df.rename(columns={pos_col:"POS"}, inplace=True)
        if settingsClass.use_z == True:
            raw_df["P"] = stats.norm.sf(abs(raw_df["Z"]))*2 # Use 2 tailed Z test to calculate p-value
            pval_col = "P"
        print("OK")

        if argumentClass.snp_map != None:
            print("Mapping SNP IDs to genomic coordinates...",end="",flush=True)
            raw_df = snpToLoc(raw_df, argumentClass.snp_map)
        if argumentClass.ref_col != None:
            raw_df = raw_df.rename(columns={argumentClass.ref_col:"REF"})
        if argumentClass.alt_col != None:
            raw_df = raw_df.rename(columns={argumentClass.alt_col:"ALT"})
    print("OK")

    print("Filtering SNPs for analysis...",end="", flush=True)
    raw_df.columns = raw_df.columns.str.replace('[#,@,&]', '',regex=True)

    pd.options.mode.chained_assignment = None # Because this chained assignment is ok, we will suppress this warning
    gwas_df = raw_df[(raw_df["CHR"]!="X") & (raw_df["CHR"]!="Y") & (raw_df["CHR"]!="M") & (raw_df["CHR"]!=23) & (raw_df["CHR"]!="23")] # Autosomes only
    gwas_df['CHR'] = gwas_df['CHR'].astype(int)
    gwas_df['POS'] = gwas_df['POS'].astype(int)
    gwas_df = gwas_df.sort_values(axis=0, by=["CHR","POS"], ascending=True)
    gwas_df[argumentClass.p_col] = gwas_df[argumentClass.p_col].replace(to_replace=list([NaN,"-"]),value=1)
    gwas_df[argumentClass.p_col] = gwas_df[argumentClass.p_col].astype(float)
    gwas_df = gwas_df[gwas_df[argumentClass.p_col]<argumentClass.pval] # p-val threshold

    if argumentClass.maf_col in gwas_df.columns: # MAF threshold
        gwas_df = gwas_df[gwas_df[argumentClass.maf_col]>argumentClass.maf]
    gwas_df = gwas_df[((gwas_df["REF"].isin(["A","T","C","G"])) & (gwas_df["ALT"].isin(["A","T","C","G"])))] # Biallelic only
    last_col = list(["ID","QUAL","FILTER","INFO","REF","ALT"])
    for c in last_col:
        if c not in gwas_df.columns:
            gwas_df[c] = "."
        else:
            gwas_df[c] = gwas_df[c].fillna(".")

    if len(gwas_df)<1:
        raise RuntimeError(f"No genome-wide significant SNPs found. SNPs must be autosomal, biallelic, and have a p-value greater than the threshold (p>{argumentClass.pval}) to be passed on for annotation. You may change this p-value threshold using the --pval option. \nExiting.")
    # if len(gwas_df)<1: # Check that we have identified genome-wide significant SNPs which meet the QC conditions
    #     print(f"No genome-wide significant SNPs found. SNPs must be autosomal, biallelic, and have a p-value greater than the threshold (p>{pval}) to be passed on for annotation. You may change this p-value threshold using the --pval option. \nExiting.")
    #     exit()

    snps_out = f"{argumentClass.output}.snps"
    gwas_df[["CHR","POS","ID","REF","ALT","QUAL","FILTER","INFO"]].to_csv(snps_out, sep="\t", index=False, header=False)

    print("OK")
    return snps_out

###############################################################################
############################# Annotation with VEP #############################
###############################################################################
def check_build(vep_file, argumentClass):
    '''
    Confirm that the VEP annotation has been created using the correct genome
    build. Failure to do so will generate incorrect results. Also flag
    if VEP API and cache version don't match.
    '''
    with open(vep_file, "r") as f:
        for line in f:
            if line.startswith('#'):
                if re.search("assembly version", line):
                    if ("GRCh"+argumentClass.build) in line:
                        return
                    else:
                        raise RuntimeError(f"VEP annotation is not using GRCh{argumentClass.build} genome reference build. Check that you are providing the correct --ref_build, --cache_dir, and --cache_version. \nExiting.")
                if re.search("API version", line):
                    if (argumentClass.cache_v) in line:
                        return
                    else:
                        print(f"\nWARNING: VEP API version does not match cache version ({argumentClass.cache_v}). VEP recommends using the same cache and API version for best results.")
            else:
                return
    return "pass"

def vep_annotate(snps_out, argumentClass):
    '''
    Use Ensembl Variant Effect Predictor (VEP) to identify features which are
    within a specified distance of the genome-wide significant SNPs. Allowed
    features are specified with --biotype and only results within those
    categories will be kept.
    '''
    print("Annotating SNPs...", end="", flush=True)
    vep_out = f"{argumentClass.output}.annotation"
    vep_summ = f"{argumentClass.output}.annotation_summary.html"

    # Check whether VEP annotation already exists and remove it with warning if True
    if os.path.exists(vep_out) == True:
        print(f'\nWARNING: {vep_out} already exists and will be removed/overwritten.')
        os.remove(vep_out)
        if os.path.exists(vep_summ) == True:
            os.remove(vep_summ)

    # Now run VEP and store output in .annotation file
    ann_names = ["VAR_PROVIDED", "LOCATION", "ALLELE", "GENE", "FEATURE", "FEATURE_TYPE", "CONSEQUENCE", "cDNA_POS", "CDS_POS", "PROT_POS", "AMINO_ACIDS", "CODONS", "EXISTING_VARIATION", "EXTRA"]
    if argumentClass.biotypes == "regulatory":
        cmd_annotate = ['vep', '--offline', '--cache', '--dir_cache', argumentClass.cache_dir, '--assembly', ('GRCh'+argumentClass.build), '-i', snps_out, '-o', vep_out, '--cache_version', argumentClass.cache_v, '--distance', str(argumentClass.dist), '--biotype','--force_overwrite','--regulatory']
    else:
        cmd_annotate = ['vep', '--offline', '--cache', '--dir_cache', argumentClass.cache_dir, '--assembly', ('GRCh'+argumentClass.build), '-i', snps_out, '-o', vep_out, '--nearest', 'gene', '--cache_version', argumentClass.cache_v, '--distance', str(argumentClass.dist), '--biotype','--force_overwrite']
    run_annotate = subprocess.run(cmd_annotate, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Confirm that file was created
    try:
        annotation_df = pd.read_csv(vep_out, sep="\t", skipinitialspace=True, index_col=False, names=ann_names, usecols=["LOCATION","ALLELE","GENE","FEATURE","FEATURE_TYPE","EXTRA"], comment='#')
    except FileNotFoundError:
        cmd_msg = " ".join(cmd_annotate)
        print(f'\nFileNotFoundError: VEP annotation output file not found.\n\
        This often occurs when there are issues with the VEP cache or version compatibility.\n\
        To troubleshoot, try running the command: {cmd_msg} \n\
        Exiting.')
        exit()

    # Confirm that file used the right genome reference assembly version
    check_build(vep_out, argumentClass)

    # Confirm that the file is not empty
    if len(annotation_df)<1:
        raise RuntimeError(f'VEP annotation failed. This often occurs due to issues with the VEP cache, version compatibility, or column headers. \nTo troubleshoot, try running the command: {cmd_msg} \nExiting.')

    # Remove duplicate genes and limit to only genes in the desired biotype
    annotation_df[["X","BIOTYPE"]] = annotation_df["EXTRA"].str.split("BIOTYPE=", expand=True)
    if argumentClass.biotypes == "protein_coding" or argumentClass.biotypes == "protein_all" or argumentClass.biotypes == "all":
        annotation_df[["BIOTYPE","NEAREST"]] = annotation_df["BIOTYPE"].str.split(";NEAREST=", expand=True)
    if argumentClass.biotypes == "protein_coding": # Protein coding is only allowed biotype
        allowed_biotypes=list(["protein_coding"])
        annotation_df = annotation_df[annotation_df["BIOTYPE"].isin(allowed_biotypes)]
        dup_drop = "NEAREST"
    elif argumentClass.biotypes == "protein_all": # All protein-related biotypes
        allowed_biotypes=list(["protein_coding","IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_M_gene","IG_V_gene","IG_Z_gene","nonsense_mediated_decay","nontranslating_CDS","non_stop_decay","polymorphic_pseudogene","TR_C_gene","TR_D_gene","TR_J_gene"])
        annotation_df = annotation_df[annotation_df["BIOTYPE"].isin(allowed_biotypes)]
        dup_drop = "NEAREST"
    elif argumentClass.biotypes == "all": # Do not restrict biotypes
        dup_drop = "NEAREST"
    elif argumentClass.biotypes == "regulatory": # Use only regulatory features
        annotation_df = annotation_df[annotation_df["FEATURE_TYPE"]=="RegulatoryFeature"]
        dup_drop = "FEATURE"
    annotation_df = annotation_df.drop_duplicates(subset=dup_drop).drop(columns=["X","EXTRA"])

    # Check that there are annotations left over after pruning biotypes
    if len(annotation_df)<1:
        raise RuntimeError("No annotations of compatible biotypes found. You can change the allowed biotypes using the --biotypes option. Exiting.")

    annotation_out = f"{argumentClass.output}.features"
    annotation_df.to_csv(annotation_out, sep="\t", index=False, header=True)
    print("OK")
    return annotation_out

###############################################################################
################### Gene and location finding with BioMart ####################
###############################################################################
def biomart_locate(annotation_out, argumentClass):
    '''
    Using BioMart wget queries, retrieve the genomic location (CHR, START, END)
    for each element in the annotation list. These will allow us to look for
    intersections with the known element set later.
    '''
    print("Locating elements with BioMart...",end="", flush=True)
    # Set proper output header depending on which biotypes we are looking at
    if argumentClass.biotypes == "regulatory":
        feature_col = "FEATURE"
        biomart_header = ["CHR","START","END","FEATURE_TYPE","FEATURE_ID"]
    else:
        feature_col = "NEAREST"
        biomart_header = ["ENSEMBL_ID","START","END","CHR","GENE_NAME","STRAND","HGNC_SYMBOL"]

    # Grab the location (CHR, POS) of the genes identified through vep_annotate()
    annotation_ids = pd.read_csv(annotation_out, sep="\t", header=0, usecols=[feature_col])
    biomart_out = f"{argumentClass.output}.biomart"
    biomart_df = pd.DataFrame(columns=biomart_header)
    biomart_tmp = f"{biomart_out}.tmp"
    ids_count = len(annotation_ids[feature_col])

    # Check if BioMart output file already exists and remove with warning if so
    if os.path.exists(biomart_out) == True:
        print(f'\nWARNING: {biomart_out} already exists and will be removed/overwritten.')
        os.remove(biomart_out)
        if os.path.exists(biomart_tmp) == True:
            os.remove(biomart_tmp)

    if argumentClass.build == "38":
        ensembl_build = "ensembl.org"
    else:
        ensembl_build = "grch37.ensembl.org"

    # Biomart doesn't like to accept lists longer than 400, so we will break it up if this is the case
    if ids_count<=400:
        genes = list(annotation_ids[feature_col].fillna(""))
        ensembl_ids = ",".join(genes)
        if argumentClass.biotypes == "regulatory":
            cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_regulatory_feature\" interface = \"default\" ><Filter name = \"regulatory_stable_id\" value = \"{ensembl_ids}\"/><Attribute name = \"chromosome_name\" /><Attribute name = \"chromosome_start\" /><Attribute name = \"chromosome_end\" /><Attribute name = \"feature_type_name\" /><Attribute name = \"regulatory_stable_id\" /></Dataset></Query>\'"
        else:
            cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"ensembl_gene_id\" value = \"{ensembl_ids} \"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"strand\" /><Attribute name = \"hgnc_symbol\" /></Dataset></Query>\'"
        os.system(cmd_biomart)
        tmp_df = pd.read_csv(biomart_tmp, sep="\t", names=biomart_header)
        biomart_df = pd.concat([biomart_df, tmp_df], axis=0, ignore_index=True)

    else:
        range_list = list(range(400,ids_count,400))
        range_list.insert(0,0)
        if range_list[-1] != ids_count:
            range_list.append(ids_count)
        for i in range(len(range_list)-1):
            start = range_list[i]
            end = range_list[i+1]
            genes = list(annotation_ids[feature_col].fillna(""))[start:end]
            ensembl_ids = ",".join(genes)
            if argumentClass.biotypes == "regulatory":
                cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_regulatory_feature\" interface = \"default\" ><Filter name = \"regulatory_stable_id\" value = \"{ensembl_ids}\"/><Attribute name = \"chromosome_name\" /><Attribute name = \"chromosome_start\" /><Attribute name = \"chromosome_end\" /><Attribute name = \"feature_type_name\" /><Attribute name = \"regulatory_stable_id\" /></Dataset></Query>\'"
            else:
                cmd_biomart = f"wget -O {biomart_tmp} -q \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"ensembl_gene_id\" value = \"{ensembl_ids} \"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"strand\" /><Attribute name = \"hgnc_symbol\" /></Dataset></Query>\'"
            os.system(cmd_biomart)
            tmp_df = pd.read_csv(biomart_tmp, sep="\t", names=biomart_header)
            biomart_df = pd.concat([biomart_df, tmp_df], axis=0, ignore_index=True)

    # Confirm that BioMart location finding succeeded and returned results
    if len(biomart_df) < 1:
        cmd_msg = " ".join(cmd_biomart)
        raise RuntimeError(f"\nRuntimeError: BioMart location finding failed.\n\
        Please check that \'wget\' is installed in your environment and that your reference file is unzipped.\n\
        You can attempt to troubleshoot BioMart using the command {cmd_msg}.\n\
        If this error persists, please open an issue on GitHub. \nExiting.")
        exit()
    biomart_df.to_csv(biomart_out, sep="\t", header=True, index=False)

    print("OK")
    print(f"[biomart] Total number of annotations: {ids_count}")

    # Reformat to BED file format for import into bedtools and EOI enrichment testing
    locations_out = f"{argumentClass.output}.locations.bed"
    print(f"Writing output to {locations_out}...", end="", flush=True)
    locations_df = pd.DataFrame({})
    if argumentClass.build == "38":
        biomart_df["CHR"] = "chr" + biomart_df["CHR"].astype(str)
    locations_df[1], locations_df[2], locations_df[3] = biomart_df["CHR"], biomart_df["START"], biomart_df["END"]
    locations_df.to_csv(locations_out, sep="\t", header=False, index=False)
    print("OK")
    return locations_out

###############################################################################
##################### Simulation and enrichment analysis ######################
###############################################################################
def sim_prep(locations_out):
    '''
    Load element set for enrichment analysis and compute statistics on them.
    The data about element set size and lengths will be used to create a
    matched element set when performing simultions.
    '''
    print("Loading element set...",end="", flush=True)
    df_elemSet = pd.read_csv(locations_out, skipinitialspace=True, sep="\t", header=None, names=["CHR","START","END"], usecols=[0,1,2])
    print("OK")

    # Get details about the element set that we can use to make matched set(s)
    sSize = len(df_elemSet.index)
    eDiff = list()
    df_elemSet["LENGTH"] = abs(df_elemSet["END"]-df_elemSet["START"])
    df_elemSet["BIN"] = pd.qcut(df_elemSet["LENGTH"],q=10,duplicates="drop")
    bin_len = df_elemSet.groupby(["BIN"], observed=False).mean()["LENGTH"].reset_index()
    bin_size = df_elemSet.groupby(["BIN"], observed=False).count()["LENGTH"].reset_index() # All columns have the same value, which is a count of how many elements are in the bin
    bin_bp = df_elemSet["LENGTH"].sum()
    lens = list()
    sizes = list()
    for i in range(len(bin_len)):
        lens.append(np.floor(bin_len["LENGTH"][i]))
        sizes.append(bin_size["LENGTH"][i])
    print("[elements] Average element lengths: ", lens)
    print("[elements] Number of elements in each bin: ", sizes)
    print(f"[elements] Element set size: {bin_bp} bp")
    return bin_len, bin_size, bin_bp, sSize

###############################################################################
###############################################################################
###############################################################################
# def sim_matchFile():
#     simulation_outPath = f"{output}.simulation.tmp"
#     create deciles
#     for d in deciles:
#         cmd_match = ["bedtools"]
#     with open(simulation_outPath, "w") as outfile:
#     for i in range(0,len(bin_len)):
#         L = str(bin_len["LENGTH"][i])
#         S = str(bin_size["LENGTH"][i])
#
#         if use_seed == True: # Only for unittest
#             cmd_rand = ["bedtools","random","-l",L,"-n",S,"-g",argumentClass.ref,"-seed","1001"]
#         else:
#             cmd_rand = ["bedtools","random","-l",L,"-n",S,"-g",argumentClass.ref]
#         # print(cmd_rand)
#         run_rand = subprocess.run(cmd_rand, stdout=outfile, stderr=subprocess.PIPE)
# outfile.close()
#
# if os.stat(simulation_outPath).st_size == 0:
#     cmd_msg = " ".join(cmd_rand)
#     raise RuntimeError(f"RuntimeError: bedtools simulation failed. Please check your reference file.\n\
#     You may also troubleshoot using this command: {cmd_msg}\n\
#     Exiting.")
#     exit()
# return simulation_outPath
###############################################################################
###############################################################################
###############################################################################

def simulate(bin_len, bin_size, argumentClass, use_seed):
    '''
    In order to test for enrichment, use BEDTools to simulate a set of random
    gene regions using the reference genome annotation. These will be used to
    create a background distribution that we can perform a statistical test
    on to determine enrichment.
    '''
    simulation_outPath = f"{argumentClass.output}.simulation.tmp"

    # options = match_on(whatever_here)

    with open(simulation_outPath, "w") as outfile:
        for i in range(0,len(bin_len)):
            L = str(bin_len["LENGTH"][i])
            S = str(bin_size["LENGTH"][i])

            if use_seed == True: # Only for unittest
                cmd_rand = ["bedtools","random","-l",L,"-n",S,"-g",argumentClass.ref,"-seed","1001"]
            else:
                cmd_rand = ["bedtools","random","-l",L,"-n",S,"-g",argumentClass.ref]
            # print(cmd_rand)
            run_rand = subprocess.run(cmd_rand, stdout=outfile, stderr=subprocess.PIPE)
    outfile.close()

    if os.stat(simulation_outPath).st_size == 0:
        cmd_msg = " ".join(cmd_rand)
        raise RuntimeError(f"RuntimeError: bedtools simulation failed. Please check your reference file.\n\
        You may also troubleshoot using this command: {cmd_msg}\n\
        Exiting.")
        exit()
    return simulation_outPath

def intersect(argumentClass, fileB, eBP):
    '''
    Find intersections between element set and elements of interest.
    Our enrichment analysis uses the statistic of intersections per base pair.
    '''

    intersect_outPath = f"{argumentClass.output}.intersect.tmp"
    cmd_int = ["bedtools","intersect","-c","-a",str(fileB),"-b",str(argumentClass.eoi)]
    with open(intersect_outPath, "w") as outfile:
        run_int=subprocess.Popen(cmd_int,stdin=subprocess.PIPE,stdout=outfile)
        stdout, stderr = run_int.communicate()

    if os.stat(intersect_outPath).st_size == 0:
        print(f"RuntimeError: bedtools intersect failed. Please check your reference file.\n\
        You may also troubleshoot using this command: {cmd_int}\n\
        Exiting.")
        exit()

    with open(intersect_outPath,"r") as f:
        int_count = 0
        for line in f:
            int_count += int(line.split()[-1])

    int_per_bp = int_count/eBP
    return int_per_bp

###############################################################################
#################################### MAIN #####################################
###############################################################################
def main(**kwargs):
    # Declare settings and parameters coming from options/arguments to hare.py
    hareSettings = hareclasses.SettingsContainer(kwargs["use_z"], kwargs["anno_only"], kwargs["source_neale"],
    kwargs["source_bolt"], kwargs["keep_tmp"]) # kwargs["not_human"]

    hareParameters = hareclasses.ArgumentContainer(kwargs["gwas"], kwargs["pval"], kwargs["gwas_p"],
    kwargs["maf"], kwargs["gwas_maf"], kwargs["gwas_ref"], kwargs["gwas_alt"], kwargs["snp_map"],
    kwargs["out"], kwargs["cache_dir"], kwargs["cache_version"], kwargs["biotypes"], kwargs["dist"],
    kwargs["eoi"], kwargs["ref"], kwargs["ref_build"], kwargs["draws"], hareSettings)

    print("\n ----------------------------------------------------------------------")
    print("|                                                                      |")
    print("|      Welcome to the HARE genetic feature enrichment pipeline.        |")
    print("|  Please contact Olivia Smith at osmith@utexas.edu to report issues.  |")
    print("|                                                                      |")
    print(" ----------------------------------------------------------------------\n")

    print(f"-----Workflow started at {datetime.now()}.\n")
    print(f"[args] GWAS summary statistics file: \'{hareParameters.gwas}\'")
    if hareSettings.source_neale == True:
        print(f"[args] Using Neale Lab summary statistics format.")
    elif hareSettings.source_bolt == True:
        print(f"[args] Using BOLT-LMM summary statistics format.")
    if hareSettings.use_z == True:
        print("[args] --use_z is ON. Will compute p-values from two-tailed Z test.")
    print(f"[args] P-value significance threshold: {hareParameters.pval}")
    print(f"[args] Minor allele frequency (MAF) threshold: {hareParameters.maf}")
    print(f"[args] VEP cache directory: \'{hareParameters.cache_dir}\', version {hareParameters.cache_v}")
    print(f"[args] Annotation max distance: {hareParameters.dist}")
    print(f"[args] Biotypes included in analysis: {hareParameters.biotypes}")
    print(f"[args] Reference genome annotation file: \'{hareParameters.ref}\'")
    # if hareSettings.non_human == True:
    #     print(f"[args] WARNING: Using --non-human option. This is currently allowed but not explicitly supported.")
    #     print(f"[args] Reference genome build: {hareParameters.build}")
    # else:
    #     print(f"[args] Reference genome build: GRCh{hareParameters.build}")
    print(f"[args] Reference genome build: GRCh{hareParameters.build}")
    if hareSettings.anno_only == False:
        print(f"[args] Element of interest BED file: \'{hareParameters.eoi}\'")
        print(f"[args] Simulation draws to be performed: {hareParameters.draws}")
    else:
        print(f"[args] Annotate only is ON, will not perform simulation and intersection analysis.")
    print(f"[args] Results will be written to {hareParameters.output}.*")

    # Check that all the files and directories exist
    for f in [hareParameters.gwas, hareParameters.snp_map, hareParameters.cache_dir, hareParameters.eoi, hareParameters.ref]:
        findFiles(f)

    # Run annotation and feature finding
    snps_filepath = gwas_import(hareParameters, hareSettings)
    annotation_filepath = vep_annotate(snps_filepath, hareParameters)#, hareSettings)
    locations_filepath = biomart_locate(annotation_filepath, hareParameters)#, hareSettings)
    if hareSettings.anno_only == True:
        print(f"\nWorkflow completed at {datetime.now()}.-----\U0001F407\n")
        return

    # Intersection testing
    intersections = list()
    cat = list()
    len_vector, s_vector, total_bp, set_size = sim_prep(locations_filepath)

    # Generate and test simulation set
    print("Simulating matched element sets and identifying intersections...", end="", flush=True)
    for d in range(hareParameters.draws):
        # Introduce function to match with recombination dist
        # simulation_outPath = f"{output}.simulation.tmp"
        # if matchFile != None:
        #     simulation_results = sim_matchFile(PARAMS, simulation_outPath)
        # else:
        #     simulation_results = simulate(len_vector, s_vector, ref, simulation_outPath)
        simulation_results = simulate(len_vector, s_vector, hareParameters, False)
        intersections.append(intersect(hareParameters, simulation_results, total_bp))
        cat.append("simulation")
    print("OK")

    # Test element set of interest
    print("Examining test element set...",end="", flush=True)
    set_intersect = intersect(hareParameters, locations_filepath, total_bp)
    intersections.append(set_intersect)
    cat.append("test_set")
    print("OK")

    #### Generate output file ####
    print("Generating output file...",end="", flush=True)
    out_df = pd.DataFrame({"category":cat, "int_per_bp":intersections, "set_size":set_size})
    out_path = f"{hareParameters.output}.intersections"
    out_df.to_csv(out_path,sep="\t",index=False)
    print("OK")

    ### Remove temporary files ###
    if hareSettings.keep_tmp != True:
        print("Cleaning up directory...",end="", flush=True)
        os.remove(f"{hareParameters.output}.biomart.tmp")
        os.remove(simulation_results)
        os.remove(f"{hareParameters.output}.intersect.tmp")
        print("OK")

    print(f"\nWorkflow completed at {datetime.now()}.-----\n")
    return

if __name__ == "__main__":
    main()
