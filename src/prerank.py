#!/usr/in/env python
# enrichmentPrep.py
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

################################################################################

# This script creates a ranked list file of genes based on an input p-value file.
# This script is part of the HARE toolkit. These ranked score files are
# compatible with tools like GSEA and WebGestalt.
# It can also be used for peak calling by identifying genes with multiple
# positions which have genome-wide significant p-values.

# Dependencies: python => 3.0 and python packages pandas, numpy, wget (for BioMart access)

# Inputs: [INPUT]: File containing p-values associated with each position.
#                   Chromosome, position, and p-value columns are required.
#         [INCLUDE].bed (optional): BED format file containing a list of
#                   positions to filter to for analysis. Only listed positions
#                   will be included.
#         [EXCLUDE].bed (optional): BED format file containing a list of
#                   positions to remove from analysis.

# Outputs: [OUTPUT].rnk: File containing HGNC symbol and score for each annotated
#                        gene. Use --score_method to choose whether min or mean
#                        p-value from associated positions is used to score.
#          [OUTPUT].dmp: File containing peaks for which a sufficient number of
#                        positions within a given distance meet genome-wide
#                        significance threshold. Example: peaks with at least 3
#                        positions within 1,000 bp of a gene have a p-value
#                        below 1e-06.

# Example command: hare prerank --input [INPUT] --output [OUTPUT] ... [OPTIONS]

################################################################################

###############################################################################
################ Setup libraries, arguments, and dependencies #################
###############################################################################
import pandas as pd
import numpy as np
import sys
import subprocess
import argparse
import os
from shutil import which
from datetime import datetime
import hareclasses

def findFiles(filename):
    if filename == None:
        return
    if os.path.exists(filename) == False:
        raise FileNotFoundError(f'{filename} does not exist or could not be opened.')
    return

def checkHeader(argumentClass):
    print(argumentClass.vCol, argumentClass.cCol, argumentClass.pCol)
    tableHead = pd.read_csv(argumentClass.input, sep="\t", header=0, nrows=3)
    requiredCols = [argumentClass.vCol, argumentClass.cCol, argumentClass.pCol]
    for h in requiredCols:
        if h not in tableHead.columns:
            raise KeyError(f"Necessary column, {h}, is not present in input file. Columns specifying chromosome, position, and p-value all required. \nSpecify different column headers using --<HEADER>_col.")
    return

###############################################################################
######################## Filter positions for analysis ########################
###############################################################################

def filterPositions(argumentClass):
    '''
    Filter positions based on p-value and any filtering (inclusion or exclusion)
    BED files. Can also constrain to specific number of top hits or chromosome.
    '''

    print("Reading in p-values file...", end="", flush=True)
    pvalDF = pd.read_csv(argumentClass.input, sep="\t", header=0)
    pvalDF.rename(columns={argumentClass.vCol:"P", argumentClass.cCol:"CHR", argumentClass.pCol:"POS"}, inplace=True)
    print("OK")

    print("Filtering positions...", end="", flush=True)

    # Filter to only positions specified in the provided BED file
    if argumentClass.fFile != None:
        fDF = pd.read_csv(argumentClass.fFile, sep="\t", header=None)
        fDF.rename(columns={0: "CHR", 1:"POS", 2:"END"}, inplace=True)
        pvalDF = pd.merge(pvalDF, fDF, on=["CHR", "POS"], how="inner")

    # Remove any positions that were specified in the provided BED file
    if argumentClass.xFile != None:

        xDF = pd.read_csv(argumentClass.xFile, sep="\t", header=None)
        xDF.rename(columns={0: "CHR", 1:"POS", 2:"END"}, inplace=True)

        if argumentClass.fFile == None:
            pvalDF = pd.merge(pvalDF, xDF, on=["CHR", "POS"], how="outer", indicator=True).query('_merge == "left_only"').drop(columns='_merge')
        else:
            pvalDF = pd.merge(pvalDF, xDF, on=["CHR", "POS", "END"], how="outer", indicator=True).query('_merge == "left_only"').drop(columns='_merge')

            # There should be no overlap between included and excluded positions
            if pd.merge(xDF, fDF, on=["CHR","POS","END"], how="inner").empty != True:
                raise ValueError("Files designating positions to include and exclude have overlapping positions.")

    pvalDF["END"] = pvalDF["POS"].astype(int) + 1

    if argumentClass.chrReg != None:
        pvalDF = pvalDF[pvalDF["CHR"]==argumentClass.chrReg] # Restrict to given chromosome

    # Use all positions if no file given, otherwise get N most significant values
    if argumentClass.nTop == None:
        pass
    else:
        pvalDF = pvalDF.nsmallest(argumentClass.nTop, "P")

    if argumentClass.pThresh == 1:
        pass
    else:
        pvalDF = pvalDF[pvalDF["P"].astype(float)<=argumentClass.pThresh]

    if len(pvalDF) == 0:
        raise RuntimeError("No positions left after filtering. Check any inclusion/exclusion filters, chromosome restrictions, and p-value thresholds.")

    pvalDF["POS"] = pvalDF["POS"].astype(int)
    pvalDF.sort_values(by = ["CHR","POS"], ignore_index=True)
    pvalDF = pvalDF.reset_index(drop=True)

    # pvalDF.to_csv(f"{argumentClass.output}.pvalDF.txt", sep="\t", header=True, index=False)

    print("OK")
    return pvalDF

###############################################################################
########################### Annotate with BioMart #############################
###############################################################################

def annotate(pvalDF, argumentClass):
    '''
    Annotate genes using BioMart which are within a buffer distance of upstream/
    downstream. Can use either Ensembl GRCh38 or GRCh37. Gets HGNC symbols which
    will be used for the ranking.
    '''
    print("Annotating genes...")

    biomartTMP = f"{argumentClass.output}.biomart.tmp"
    outPath = f"{argumentClass.output}.biomart"

    if os.path.exists(outPath) == True:
        print(f'\nWARNING: {outPath} already exists and will be removed/overwritten.')
        os.remove(outPath)
        if os.path.exists(biomartTMP) == True:
            os.remove(biomartTMP)

    headList = ["ENSEMBL_ID", "CHR", "GENE_START", "GENE_END", "GENE_SYMBOL", "HGNC_SYMBOL", "BIOTYPE"]
    biomartDF = pd.DataFrame(columns=headList)

    chrList = list(pvalDF["CHR"])
    posList = list(pvalDF["POS"])
    # scoreList = list(logDF["logP"])
    chrPosEnd = list()
    for i in range(len(posList)):
        chrPosEnd.append(f"{chrList[i]}:{posList[i]-argumentClass.buffer}:{posList[i]+argumentClass.buffer}")

    if argumentClass.build == "37":
        ensembl_build = "grch37.ensembl.org"
    elif argumentClass.build == "38":
        ensembl_build = "ensembl.org"
    else:
        raise KeyError("Ensembl build must either be 37 or 38.")

    if len(chrPosEnd)<=250:
        inList = ",".join(str(p) for p in chrPosEnd)
        cmd_biomart = f"wget -O {outPath} \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?> <!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"chromosomal_region\" value = \"{inList}\"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"hgnc_symbol\" /><Attribute name = \"gene_biotype\" /></Dataset></Query>\'"
        os.system(cmd_biomart)
        biomartDF = pd.read_csv(outPath, sep="\t", names=headList)
    else:
        range_list = list(range(250,len(chrPosEnd),250))
        range_list.insert(0,0)
        if range_list[-1] != len(chrPosEnd):
            range_list.append(len(chrPosEnd))
        for i in range(len(range_list)-1):
            r_start = range_list[i]
            r_end = range_list[i+1]
            inList = ",".join(str(p) for p in chrPosEnd[r_start:r_end])
            cmd_biomart = f"wget -O {biomartTMP} \'https://{ensembl_build}/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?> <!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"chromosomal_region\" value = \"{inList}\"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"chromosome_name\" /><Attribute name = \"start_position\" /><Attribute name = \"end_position\" /><Attribute name = \"external_gene_name\" /><Attribute name = \"hgnc_symbol\" /><Attribute name = \"gene_biotype\" /></Dataset></Query>\'"
            try:
                os.system(cmd_biomart)
            except:
                raise RuntimeError(f"BioMart annotation failed. Troubleshoot with the command:\n\n{cmd_biomart}\n")
            time.sleep(0.5) # Don't want to overwhelm the BioMart system because it will stop allowing connections
            try:
                tmpDF = pd.read_csv(biomartTMP, sep="\t", names=headList)
            except:
                raise RuntimeError(f"BioMart annotation failed. Troubleshoot with the command:\n\n{cmd_biomart}\n")
            biomartDF = pd.concat([biomartDF, tmpDF], axis=0, ignore_index=True)
            # print(biomartDF)
    biomartDF.dropna(subset=["HGNC_SYMBOL"], inplace=True) # HGNC symbols required for these tools
    if argumentClass.biotypes == 'coding':
        biomartDF = biomartDF[biomartDF["BIOTYPE"] == 'protein_coding']

    biomartDF[["CHR", "GENE_START", "GENE_END"]] = biomartDF[["CHR", "GENE_START", "GENE_END"]].apply(pd.to_numeric)
    biomartDF = biomartDF.sort_values(by=["CHR","GENE_START"], ascending=[True, True]).reset_index(drop=True)
    biomartDF.to_csv(outPath, sep="\t", header=True, index=False) # Write final version to file

    print("Annotation complete.")
    return biomartDF

###############################################################################
################### Score genes and call significant peaks ####################
###############################################################################

def call_peaks(logDF, biomartDF, argumentClass):
    '''
    Call peaks where at least N positions associated with a given gene (again
    within a buffer window up and downstream of that gene) are below the p-value
    threshold.
    '''
    print("Calling peaks...", end="", flush=True)

    biomartDF['SIG_COUNT'] = 0
    for index, row in biomartDF.iterrows():
        count = 0
        filtered_rows = logDF[(logDF['CHR'] == row['CHR']) & (logDF['P'] < argumentClass.dmpP)]
        up_bound = row['GENE_START'] - argumentClass.dmpD
        down_bound = row['GENE_END'] + argumentClass.dmpD

        for position in filtered_rows["POS"]:
            if (position >= up_bound) & (position <= down_bound):
                count = count + 1
        biomartDF.loc[index, 'SIG_COUNT'] = count

    biomartDF[(biomartDF["SIG_COUNT"]>=argumentClass.dmpThresh)].to_csv(f"{argumentClass.output}.dmp", sep="\t", header=True, index=False)
    print("OK")
    return biomartDF

def score_and_rank(pvalDF, biomartDF, argumentClass):
    '''
    Generate score for each gene using -log10(p) values. Can either use the
    minimum or average score of all associated positions for the final score
    for the gene. Genes will be ranked (placed in order) after scoring.
    '''
    print("Scoring...", end="", flush=True)

    logDF = pvalDF.sort_values(by=["CHR","POS"], ascending=[True, True]).reset_index(drop=True)
    logDF["logP"] = np.log10(logDF["P"])*-1

    if argumentClass.call_peaks == True:
        biomartDF = call_peaks(logDF, biomartDF, argumentClass)

    outGeneName = list()
    outScores = list()
    for r in range(len(biomartDF["GENE_START"])):
        start = biomartDF["GENE_START"][r]
        end = biomartDF["GENE_END"][r]
        chr = biomartDF["CHR"][r]
        for s in range(len(logDF["CHR"])):
            refChr = logDF["CHR"].iloc[s]
            refPos = logDF["POS"].iloc[s]
            refScore = logDF["logP"].iloc[s]
            if refChr != chr:
                continue
            else:
                if refPos < (start - argumentClass.buffer):
                    continue
                elif (refPos >= (start - argumentClass.buffer)) & (refPos <= (end + argumentClass.buffer)):
                    outGeneName.append(biomartDF["HGNC_SYMBOL"][r])
                    outScores.append(refScore)
                else:
                    break

    if argumentClass.scoring == 'min':
        dfRanks = pd.DataFrame({"HGNC_SYMBOL":outGeneName, "SCORE":outScores}).groupby(["HGNC_SYMBOL"]).max()
    elif argumentClass.scoring == 'mean':
        dfRanks = pd.DataFrame({"HGNC_SYMBOL":outGeneName, "SCORE":outScores}).groupby(["HGNC_SYMBOL"]).mean()

    dfRanks = dfRanks.sort_values(by=["SCORE"], ascending=False)#.reset_index(drop=True)
    if argumentClass.call_peaks == True:
        dfRanks = pd.merge(dfRanks, biomartDF[["HGNC_SYMBOL", "SIG_COUNT"]], on="HGNC_SYMBOL") # Add significant called peaks count to dataframe

    fRank = f"{argumentClass.output}.rnk"
    dfRanks.to_csv(fRank, sep="\t", header=False, index=True)

    print("OK")
    return

def main(**kwargs):
    prerankArgs = hareclasses.PrerankArgumentContainer(kwargs["input"], kwargs["output"], kwargs["ref_build"],
    kwargs["buffer"], kwargs["biotypes"], kwargs["topN"], kwargs["pval_col"], kwargs["chr_col"],
    kwargs["pos_col"], kwargs["pval"], kwargs["dmpN"], kwargs["dmpP"], kwargs["dmpD"],
    kwargs["chr"], kwargs["excl"], kwargs["incl"], kwargs["score_method"], kwargs["call_peaks"])

    print("\n ----------------------------------------------------------------------")
    print("|                                                                      |")
    print("|      Welcome to the HARE genetic feature enrichment pipeline.        |")
    print("|  Please contact Olivia Smith at osmith@utexas.edu to report issues.  |")
    print("|                                                                      |")
    print(" ----------------------------------------------------------------------\n")

    fileList = [prerankArgs.input, prerankArgs.xFile, prerankArgs.fFile]
    for f in fileList:
        findFiles(f)

    # if (os.exists("wget") == False):
    #     raise RuntimeError("wget is either not installed or not executable in this environment and is required for this workflow.")

    print(f"\nEnrichment prep workflow started at {datetime.now()}.")
    print(f"\n[args] Input file: {prerankArgs.input}")
    print(f"[args] Output stem: {prerankArgs.output}")
    if prerankArgs.nTop != None:
        print(f"[args] Filtering to {prerankArgs.nTop} lowest p-values for calculation")
    print(f"[args] Filtering to p-values below {prerankArgs.pThresh}")
    print(f"[args] Allowed biotypes: {prerankArgs.biotypes}")
    if prerankArgs.chrReg == None:
        print("[args] Running on all chromsomes")
    else:
        print(f"[args] Running on chromosome {prerankArgs.chrReg}")
    if prerankArgs.xFile != None:
        print(f"[args] Excluding positions found in {prerankArgs.xFile}")
    if prerankArgs.fFile != None:
        print(f"[args] Filtering to positions found in {prerankArgs.fFile}")
    print(f"[args] Will call peaks which have {prerankArgs.dmpThresh} positions within {prerankArgs.buffer} bp below {prerankArgs.dmpP} ")

    checkHeader(prerankArgs) # Confirm that header is correct
    filteredDF = filterPositions(prerankArgs) # Filter positions
    biomartAnnotation = annotate(filteredDF, prerankArgs) # Annotate genes to positions of interest
    score_and_rank(filteredDF, biomartAnnotation, prerankArgs) # Score and rank

    print(f"\nEnrichment prep workflow completed at {datetime.now()}")

if __name__ == "__main__":
    main()
