This directory includes example input and output files for the HARE pipeline. File descriptions are provided below.

## Test run script
`example_run.sh`: This script will use the included example files to run the software. This requires the user to have 1. a genome reference and 2. a VEP cache, though users who wish to use the default ~$HOME/.vep directory for cache can remove the `--cache_dir` flag. This script will create files with the file stem test_example.* in order to prevent accidental overwriting of the provided example files.

To run:  
`bash example_run.sh {REFERENCE_PATH} {CACHE_DIRECTORY}`

## Input examples
- `example_gwas.txt`: An example of a GWAS summary statistics file which can be used as input. It includes the necessary genomic position details (chromosome, position), reference and alternate alleles, minor allele frequencies, and P-values.
- `example_eoi.bed`: Example file with regions/elements of interest (EOIs) in BED format. You can also find the human accelerated region EOIs in the main directory (`hare.reference.assets.tar.gz`).

## Output examples
- `example_result.snps`: A list of SNPs which are considered significant based on the set p-value threshold (in this case using the default value of 1e-6) and which pass QC. This file is in VCF format and used as input for the VEP annotation.
- `example_result.annotation`: Example output from VEP annotation pipeline.
- `example_result.features`: Condensed and filtered annotation. Used as input for BioMart query.
- `example_result.biomart`: Example output from BioMart. This file contains the Ensembl gene IDs and locations of the features (plus up/downstream buffer as specified).
- `example_result.intersections`: Output of `hare intersect`. The intersections/bp of each simulation as well as the test set (input element set determined from significant SNPs). **Note that this file may not be the same for every run, as this pipeline uses random simulations to create the background distribution that is intersected against**. The test_set intersections/bp value should be the same across runs.
- `example_result.png`: Output of `hare sigtest`. A .png result showing the background distribution (grey bars) and test value (red line) of intersections/bp between the simulated or test element set and the regions of interest. This distribution is generated off of n=100 draws to keep the test run short, but we do not recommend less than n=1,000 for analytical runs. A sufficiently high n is necessary in order to have a distribution which can be reasonably modeled as Weibull distributed.
- `example_result.stats`: Output of `hare sigtest`. Calculation of statistics from hypothesis testing (most notably a p-value and adjusted p-value). Used to determine whether the intersections/bp of the test set are elevated as compared to the simulations more than would be expected by chance.
- `example_result.rnk`: Output of `hare prerank`. File which contains a gene and associated score (mean -log10(p) for SNPs which are associated with that gene). This can be used for other enrichment testing pipelines like WebGestalt and GSEA prerank.
