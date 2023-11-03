#!/usr/bin/env python
# sigtest.py
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

###############################################################################

# This script performs model fitting and hypothesis testing on results generated
# using the HARE pipeline. It then generates a PNG figure showing the value of
# the intersections/bp of the phenotype-associated element set against the
# background of simulated element sets.

# Dependencies: python => 3.0, pandas, scipy, matplotlib
# See documentation for environment.yml file and links to installation documentation for dependencies

# Inputs: [RESULTS].intersections: Filename for HARE output (can also be comma-separated list) with calculated intersections/bp between HARs for randomly generated controls and provided element set

# Outputs: [OUT_STEM].stats: Tab-separated file with the test parameters and resulting p-value
#          [OUT_STEM].png: Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set

# Example command: hare sigtest --input [INPUT] --output [OUT_STEM]

###############################################################################

###############################################################################
################ Setup libraries, arguments, and dependencies #################
###############################################################################

import pandas as pd
import os
import subprocess
from datetime import datetime
import numpy as np
import ntpath
import random
from io import StringIO
from shutil import which
import matplotlib.pyplot as plt
import sys
from sys import exit

now = datetime.now().strftime("%Y%m%d%H%M%S")

###############################################################################
################### File import and significance testing ######################
###############################################################################

#### Import intersections file ####
def read_in(inFile):
    print("Loading intersection set...", end="", flush=True)

    inters = pd.read_csv(inFile, skipinitialspace=True, sep="\t", header=0)
    pd.options.mode.chained_assignment = None # Ignore chained assignments warning
    inters["int_per_bp"].loc[inters["int_per_bp"]==0] = 1e-80 # Replace 0s for significance testing compatibility
    test = inters[inters["category"]=="test_set"]
    sims = inters[inters["category"]=="simulation"]
    set_size = test.iloc[0,2]

    print("OK")
    return sims, test, set_size

#### Perform significance testing ####
def sigtest(inFile, setDist):
    '''
    Use the provided intersection data from simulations and the test to determine
    significance (p-value). Use either the provided distribution or, if 'best' is
    specified, use a goodness of fit test to find the best fit.
    '''

    print("Performing fitting and significance testing...", end="", flush=True)
    # Find ptest in the HARE/src/ directory regardless of where we run script from
    ptest_dir = os.path.dirname(__file__)
    ptest_path = os.path.join(ptest_dir, 'ptest.R')
    cmd = ['Rscript', ptest_path, '-i', inFile, '-d', setDist] # So that we can call the Rscript properly
    outR = subprocess.check_output(cmd)
    stats = outR.decode('utf-8').split(' ')
    print("OK")

    return stats

def plot(data, vert, p, title, count, output_stem):
    '''
    Plot a histogram of the background distribution and phenotype-associated
    element set value. One plot will be created for each phenotype.
    '''
    x = data["int_per_bp"]
    # Determine number of bins using Freedman–Diaconis rule
    q25, q75 = np.percentile(x, [25, 75])
    bin_width = 2 * (q75 - q25) * len(x) ** (-1/3)
    bins = round((x.max() - x.min()) / bin_width)

    # Plot histogram and add line representing test value
    plt.figure(figsize = (10, 5))
    plt.ylabel('Count')
    plt.xlabel('Intersections/BP')
    plt.title(title)
    plt.hist(x, bins=bins, histtype='bar', align='mid', orientation='vertical', color="#BAB0AC", edgecolor="#79706E")
    plt.axvline(x = vert["int_per_bp"].iloc[0], linestyle='dashed', color = '#E15759', label = 'axvline - full height')
    if count == 0:
        png_out = f"{output_stem}.png"
    else:
        png_out = f"{output_stem}_({count}).png"
    plt.savefig(png_out)
    plt.clf()

    return

###############################################################################
#################################### MAIN #####################################
###############################################################################
def main(**kwargs):

    print("\n ----------------------------------------------------------------------")
    print("|                                                                      |")
    print("|      Welcome to the HARE genetic feature enrichment pipeline.        |")
    print("|  Please contact Olivia Smith at osmith@utexas.edu to report issues.  |")
    print("|                                                                      |")
    print(" ----------------------------------------------------------------------\n")
    print(f"-----Workflow started at {datetime.now()}.\n")

    # Check inputs
    input = kwargs['input']
    distribution = kwargs['distribution']
    output_stem = kwargs['out']
    if output_stem == None:
        output_stem = "hare"
    skip_plot = kwargs['skip_plot']

    fileList = input.split(",")
    for file in fileList:
        if os.path.exists(file) == False:
            print(f'\nFileNotFoundError: {input} does not exist or could not be opened.')
            exit()
    if distribution not in ["normal", "weibull", "beta", "gamma"]:
        print(f"\nKeyError: \'--distribution\' must be from the following options: [normal, weibull, beta, gamma]. Exiting.")
        exit()

    print(f"[args] Input intersection file(s): \'{fileList}\'")
    print(f"[args] Using {distribution.upper()} distribution")
    print(f"[args] Results will be written to {output_stem}.*\n")

    outDF = pd.DataFrame(columns=["FILENAME", "SET_SIZE", "N_SIMULATIONS", "DISTRIBUTION", "KS_STAT", "SIM_IPB", "SET_IPB", "P_DISTRIBUTION", "P_EMPIRICAL"])
    c = 0
    for f in fileList:
        simDF, testDF, s = read_in(f)
        outArr = sigtest(f, distribution)
        row = [f, s] + outArr
        outDF.loc[len(outDF)] = row
        if skip_plot == False:
            plotName = os.path.splitext(ntpath.basename(f))[0]
            plot(simDF, testDF, outArr[0], plotName, c, output_stem)
            c+=1

    print("Generating output file...",end="")
    outFile = f"{output_stem}.stats"
    outDF.to_csv(outFile,sep="\t",index=False)
    print("OK")
    print(f"\nWorkflow completed at {datetime.now()}.-----\n")

if __name__ == "__main__":
    main()
