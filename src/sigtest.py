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
from sys import argv
from scipy.stats import weibull_min

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

def getStats(inFile):
    '''
    Read in and divvy up table of intersection data.
    '''
    print("Loading intersection set...", end="", flush=True)

    inters = pd.read_csv(inFile, skipinitialspace=True, sep="\t", header=0)
    pd.options.mode.chained_assignment = None # Ignore chained assignments warning
    inters["int_per_bp"].loc[inters["int_per_bp"]==0] = 1e-80 # Replace 0s for significance testing compatibility
    test = inters[inters["category"]=="test_set"]
    test_val = test.iloc[0]["int_per_bp"]
    sims = inters[inters["category"]=="simulation"]
    set_size = test.iloc[0,2]
    n_sim = len(sims["int_per_bp"])
    mean_ipb = np.mean(sims["int_per_bp"])
    # print(test_val)
    p = len(sims["int_per_bp"][sims["int_per_bp"]>test_val])/len(sims["int_per_bp"]) # Empirical p value

    stat_out = [inFile, set_size, n_sim, mean_ipb, test_val, p]
    # print(stat_out)
    print("OK")
    return sims, test, stat_out

def weibullTest(sims, test, seed):
    '''
    Use the provided intersection data from simulations and the test to determine
    empirical and Weibull-derived significance (p-value).
    '''

    print("Performing fitting and significance testing...", end="", flush=True)
    w_sims = sims["int_per_bp"].replace(0, 1e-70)
    c, loc, x = weibull_min.fit(w_sims, floc=0) # Set floc=0 because data values are positive starting from 0
    w_dist = weibull_min(c, scale=x) # Create distribution
    w_rvs = w_dist.rvs(len(w_sims), random_state=seed) # Generate RVS

    ### For testing fit visually while debugging ###
    # fig, ax = plt.subplots(1,1)
    # ax.hist(w_rvs, bins=50, density=True, alpha=0.6, color='g')
    # x_vals = np.linspace(weibull_dist.ppf(0.001), weibull_dist.ppf(0.999), 100)
    # plt.plot(x_vals, weibull_dist.pdf(x_vals), 'r-', lw=2, label='Weibull PDF')
    # ax.hist(w_sims, density=True, bins='auto', histtype='stepfilled', alpha=0.2)
    # plt.xlabel('Intersections/BP')
    # plt.ylabel('Density')
    # plt.legend()
    # plt.show()
    #################################################

    weibull_p = 1 - w_dist.cdf(test["int_per_bp"])
    print("OK")

    return [c, x, weibull_p[0]]

# def normTest(sims, test, seed):
#     '''
#     '''
#     n_dist =
#     return [c, x, normal_p[0]]

#### Perform significance testing using R script (deprecated in v1.1.0) ####
# def sigtest(inFile, setDist):
#     '''
#     Use the provided intersection data from simulations and the test to determine
#     significance (p-value). Use either the provided distribution or, if 'best' is
#     specified, use a goodness of fit test to find the best fit.
#     '''
#
#     print("Performing fitting and significance testing...", end="", flush=True)
#     # Find ptest in the HARE/src/ directory regardless of where we run script from
#     ptest_dir = os.path.dirname(os.path.realpath(argv[0]))
#     print(ptest_dir)
#     import pathlib
#     ptest_dir = os.path.dirname(__file__)
#     ptest_path = os.path.join(ptest_dir, 'ptest.R')
#     cmd = ['Rscript', ptest_path, '-i', inFile, '-d', setDist] # So that we can call the Rscript properly
#     outR = subprocess.check_output(cmd)
#     stats = outR.decode('utf-8').split(' ')
#     print("OK")
#     return stats

def plot(data, vert, title, count, output_stem):
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

    return True

def bhCorrection(p_values):
    p_sorted = np.sort(p_values) # Order p-values
    m = len(p_values) # Get number of tests
    bh_pvalues = [None] * m # Initialize adjusted p-values

    for rank in range(m): # Perform correction
        corrected = p_sorted[rank]*m/(rank+1)
        if corrected > 1.0:
            corrected = 1.0
        bh_pvalues[rank] = corrected

    return bh_pvalues

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
    distribution = 'weibull' # kwargs['distribution']
    output_stem = kwargs['out']
    if output_stem == None:
        output_stem = "hare"
    skip_plot = kwargs['skip_plot']

    fileList = input.split(",")
    for file in fileList:
        if os.path.exists(file) == False:
            print(f'\nFileNotFoundError: {input} does not exist or could not be opened.')
            exit()
    # if distribution not in ["normal", "weibull", "beta", "gamma"]:
    #     raise KeyError("'--distribution\' must be from the following options: normal, weibull, beta, gamma")
    #     exit()

    print(f"[args] Input intersection file(s): \'{fileList}\'")
    print(f"[args] Using {distribution.upper()} distribution")
    print(f"[args] Results will be written to {output_stem}.*\n")

    outDF = pd.DataFrame(columns=["FILENAME", "SET_SIZE", "N_SIMULATIONS", "SIM_IPB", "SET_IPB", "P_EMPIRICAL", "WEIBULL_SHAPE", "WEIBULL_SCALE", "P_WEIBULL"])
    file_count = 0
    for f in fileList:
        simDF, testDF, row = getStats(f)
        weibull_results = weibullTest(simDF, testDF, None) # No seed unless testing
        outDF.loc[len(outDF)] = row + weibull_results
        if skip_plot == False:
            plotName = os.path.splitext(ntpath.basename(f))[0]
            plot(simDF, testDF, plotName, file_count, output_stem)
            file_count+=1

    outDF["ADJUSTED_P"] = bhCorrection(outDF["P_EMPIRICAL"])

    print("Generating output file...",end="")
    outFile = f"{output_stem}.stats"
    outDF.to_csv(outFile,sep="\t",index=False)
    print("OK")
    print(f"\nWorkflow completed at {datetime.now()}.-----\n")

if __name__ == "__main__":
    main()
