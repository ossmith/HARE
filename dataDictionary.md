# HARE File Descriptions and Data Dictionaries
Olivia Smith  
osmith@utexas.edu  
GitHub: @ossmith

Below is comprehensive information about each of the files created using the HARE pipeline and data dictionaries for those which require it.

## HARE.py
The script which runs the analysis pipeline. There are three major parts of the workflow: 1. annotation of the genome-wide significant SNPs and creation of the element set for testing, 2. simulation of matched element sets, and 3. computation of intersections/bp of elements with human accelerated regions (HARs). See the README for details on installation and use of this script.

## HARE_Results.R
This script runs the Weibull distribution fitting and hypothesis testing of phenotype-associated element set against that distribution. See the README for details on installation and use of this script.

## environment.yml
YAML file which contains instructions for creating a conda environment named HARE which installs all the dependencies for running the analysis.

## HAR BED Files
`harsRichard2020.GRCh37.bed` and `harsRichard2020.GRCh38.bed`
BED file which contains the human accelerated regions (HARs) discovered/annotated in various different publications. This file is sourced from the supplement of [Richard, et al., 2020](https://doi.org/10.1016/j.cell.2020.02.057). The original HAR file was created with using human genome reference GRCh37. The GRCh38 HAR BED file was created using the UCSC Genome Browser liftover tool.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome |
| CHR_START | int | Starting position of the HAR (bp) |
| CHR_END  | int | End position of the HAR (bp) |
| PUB | string | Publication which identified/annotated the HAR |

## UCSC Genome Annotations
`UCSC.GRCh37.annotation.autosomes.bed.gz` and `UCSC.GRCh37.annotation.autosomes.bed.gz` are gzipped BED files with all of the gene annotations for human autosomes in GRCh37 and GRCh38. Each autosome was downloaded separately from [UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19) and then assembled into the single autosome file. Unzip this file for viewing or for use in the pipeline with the following command:
```
gunzip UCSC.GRCh37.autosomes.bed.gz
```
Descriptions are taken from the UCSC Genome Browser at time of publication.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome |
| CHR_START | int | Starting position of the feature (bp) |
| CHR_END  | int | End position of the feature (bp) |
| NAME | string | Feature name |
| SCORE | int | A score between 0 and 1000. See [UCSC Genome Browser](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format documentation for more details. |
| STRAND | string | Defines the strand. Either "." (=no strand) or "+" or "-". |
| THICK_START | int | The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, THICK_START and THICK_END are usually set to the CHR_START position. |
| THICK_END | int | The ending position at which the feature is drawn thickly (for example the stop codon in gene displays). |
| ITEM_RGB | string | An RGB value of the form R,G,B (e.g. 255,0,0). See [UCSC Genome Browser](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format documentation for more details. |
| BLOCK_COUNT | int | The number of blocks (exons) in the BED line. |
| BLOCK_SIZE | int | A comma-separated list of the block sizes. The number of items in this list should correspond to BLOCK_COUNT. |
| BLOCK_START | int | A comma-separated list of block starts. All of the BLOCK_START positions should be calculated relative to CHR_START. The number of items in this list should correspond to BLOCK_COUNT. |

## exampleGWAS.tsv
An example GWAS summary statistics file in the 'BOLT-LMM' format as used in this publication (see `--source_bolt` option); see the data dictionary for this filetype in the [BOLT-LMM v2.4 User Manual](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html). All summary statistics files must, at minimum, contain a "CHR", "POS", "P", and "MAF" column (alternative p-value and MAF column names can be specified, see `--gwas_p` and `--gwas_maf` options). HARE also accepts the Neale Lab UKB GWAS summary statistics (or files in that format, see the `--source_neale` option). You can see the dictionary for the Neale Lab GWAS analyses in [their manifest](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291).

## example.snps
A list of SNPs with genome-wide association to the phenotype in VCF file format. See [SAMtools documentation](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for more details on this file format. Details for the following data dictionary is taken from this documentation.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome associated with variant |
| POS | int | Starting position of variant |
| ID  | string | Semicolon-separated list of unique identifiers where available. If this is a dbSNP variant it is encouraged to use the rs number(s). If not provided, a value of '.' will be used. |
| REF | string | Reference base(s). Each base must be one of A,C,G,T,N (case insensitive). HARE only uses biallelic SNPs, though this is not a constraint of the VCF format. |
| ALT | string | Alternate base(s). Each base must be one of A,C,G,T,N (case insensitive). HARE only uses biallelic SNPs, though this is not a constraint of the VCF format. |
| QUAL | string | Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT is wrong). If ALT is ‘.’ (no variant) then this is −10log10 prob(variant), and if ALT is not ‘.’ this is −10log10 prob(no variant). If not provided, a value of '.' will be used. |
| FILTER | string | PASS if this position has passed all filters, i.e., a call is made at this position. Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. “q10;s50” might indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total number of samples. ‘0’ is reserved and should not be used as a filter String. If filters have not been applied, then this field should be set to the missing value. (String, no whitespace or semicolons permitted). If not provided, a value of '.' will be used. |
| INFO | string | Additional information (no whitespace, semicolons, or equals-signs permitted; commas are permitted only as delimiters for lists of values). If not provided, a value of '.' will be used. See the SAMtools VCF file format documentation for additional information, including reserved keys. |

## example.annotation
Output from the Ensembl Variant Effect Predictor command line tool. Comments are preceded by `#`, including column headers. For the file's data dictionary, see the ['default VEP output documentation'](https://uswest.ensembl.org/info/docs/tools/vep/vep_formats.html#defaultout).

## example.biomart
Output from the BioMart location finding for the elements annotated by VEP (`[OUT].annotation` file). Headers are included in this file.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| ENSEMBL_ID | string | EnsemblID for the feature |
| START | int | First position of feature |
| END  | int | Last position of feature |
| CHR | int | Chromosome the feature is located on |
| GENE_NAME | string | Gene name associated with the feature in Ensembl |
| STRAND | integer | Strand which the feature is located on. `1` is for forward and `-1` is for reverse. |

## example.locations.bed
BED file which contains only the locations of the elements annotated via VEP. This file is used as an input to the HAR intersection steps of the pipeline.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome |
| START | int | Starting position of the feature (bp) |
| END  | int | End position of the feature (bp) |

## example.intersections
BED file which contains only the locations of the elements annotated via VEP. This file is used as an input to the HAR intersection steps of the pipeline.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| category | string | Category for the calculation which specifies whether the element set is either a `simulation` or the `test_set` (phenotype-associated). |
| int_per_bp | float | Intersections per base pair computed across the entire element set. |
| set_size  | int | Number of elements present in the element set. This number should be the same across all simulations associated with a given phenotype-associated element set. |

## exampleResults.tsv
Results file which contains information about the run parameters, model fitting, and hypothesis testing (including the p-value) of the results from the main HARE script.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| file | string | Category for the calculation which specifies whether the element set is either "simulation" or "test_set" (phenotype-associated). |
| set_size | int | Number of elements present in the element set. This number applies both to the phenotype-associated and simulation element sets. |
| n_simulations  | int | Number of simulations used to generate the background distribution. |
| distribution | string | Distribution type used for fitting data for hypothesis testing. HARE_Results.R uses the Weibull distribution. |
| KS_stat | float | Kolmogorov-Smirnov (K-S) test statistic quantifying goodness-of-fit of the simulation dataset to the Weibull distribution. |
| mean_int_per_bp | float | The mean intersections/bp across all the simulations. |
| set_int_per_bp | float | The intersections/bp in the phenotype-associated element set. |
| pvalue | float | P-value from Weibull test. |
| frac_less | float | Fraction of the simulations which had lower intersections/bp than the phenotype-associated element set. Value will be between 0 and 1. |
| frac_more | float | Fraction of the simulations which had higher intersections/bp than the phenotype-associated element set. Value will be between 0 and 1. |
