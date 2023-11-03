# HARE Genetic Feature Enrichment Analysis Pipeline
Author: Olivia Smith, osmith@utexas.edu, GitHub @ossmith  
Publication Corresponding Author: Vagheesh Narasimhan, vagheesh@austin.utexas.edu, GitHub @vagheesh

This script identifies the genes and genomic features associated with autosomal genome-wide significant SNPs based on GWAS summary statistics for a given phenotype. It then tests for elevated intersections between the phenotype-associated features and the elements of interest (e.g. HARs) based on a simulated background distribution.

## Repository Structure
```
HARE
|
│   README.md    
|   dataDictionary.md
|   environment.yml
|   pyproject.toml
|   setup.cfg
|   setup.py
|   hare.reference.assets.tar.gz
|
|───src
│   │   __init__.py
│   │   hare.py
│   │   intersect.py
│   │   sigtest.py
│   │   ptest.R
│
└───example_files
|   |   exampleGWAS.tsv
|   |   example.snps
|   |   example.annotation
|   |   example.biomart
|   |   example.locations.bed
|   |   example.intersections
|   |   exampleResults.tsv
```

## Installation
### Dependencies
The following are required dependencies for HARE. Hyperlinks will direct you to the installation documentation for these dependencies. You may either install independently or through creation of a conda environment using the provided `environment.yml` file (see below).

- Linux or Unix OS
- Python => 3.0
- [VEP command line tool](https://ensembl.org/info/docs/tools/vep/script/vep_download.html) with [human reference genome cache](https://ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)
- [perl Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)
- [HTSlib](http://www.htslib.org/download/)
- [BEDtools v2.30.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- wget (Note that wget is installed by default on Linux but not Unix (Mac) operating systems, so you may need to install it separately, e.g. `brew install wget`.)
- [R => 4.1](https://www.r-project.org/) (required for significance testing)

### Install with conda and pypi
Start by cloning this repo:
```
git clone https://github.com/ossmith/HARE.git
cd HARE
```

You can install the dependencies as a conda environment using the provided `environment.yml` file. By default it will create an environment named 'hare-env'. From the HARE directory:
```
conda env create -f environment.yml
conda activate hare-env
```

You will then need to update the perl Compress::Raw::Zlib module:
```
cpan Compress::Raw::Zlib
```

And then build HARE from the included setup files with:
```
pip3 install .
```

You will then need to download a VEP genome cache if you have not already done so. Find the right cache and install using the instructions [here](https://ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). Note that the first time you run the tool after downloading the cache, it will need to do some configuration that will slow down the pipeline. This will only take place the first time you use a given cache.

Finally, you can test the install using the provided example files and reference assets. The command is:
```
tar xvf hare.reference.assets.tar.gz
gunzip reference_assets/UCSC.GRCh37.autosomes.annotation.bed.gz
hare intersect -g example_files/exampleGWAS.tsv -e reference_assets/harsRichard2020.GRCh37.bed -r reference_assets/UCSC.GRCh37.autosomes.annotation.bed -n 100 -o hare.test
hare sigtest -i hare.test.intersections
```

## Commands
### Running annotation, simulation, and intersection
The primary analysis in this repository is performed with the Python `hare intersect` function. It performs annotation, simulation, and intersection for the provided GWAS and elements of interest (for intersecting against). The most basic command is:
```
hare intersect --gwas [GWAS] --eoi [EOI] --ref [REF] --out [OUT_STEM] ... [OPTIONS]
```  
Various options are explained below.

### Running significance testing
The significance testing is done using ` sigtest`. This function is separated from the main analysis because it allows for multiple HARE output files from different GWASs to be analyzed together to create a condensed TSV. The command is:
```
hare sigtest --i [INPUT] --o [OUT_STEM] ... [OPTIONS]
```
where `--input` is the `[FILE].intersections` file (or comma-separated list of these files) created by `hare intersect`.

## File Inputs and Outputs
### intersect
Inputs  
- `[GWAS]`: Tab-separated GWAS summary statistics file which contains information on the chromosome, position, MAF, and p-value for each SNP
- `[EOI]`: BED file with elements of interest to intersect against (EOIs). The Richard, et al., 2020 HARs BED file can be found in `hare.reference.assets.tar.gz`
- `[REF]`: BED file with gene annotation of human genome reference (e.g. UCSC's hg19 .bed). UCSC genome annotations for GRCh37 and GRCh38 can be found in `hare.reference.assets.tar.gz`

Outputs  
- `[OUT_STEM].snps`: List of variants which pass QC conditions (biallelic, MAF threshold, p-value threshold) and will be used for annotation and analysis
- `[OUT_STEM].annotation`: Raw VEP annotation results
- `[OUT_STEM].annotation_summary.html`: VEP run log produced by command line run
- `[OUT_STEM].biomart`: Output from Biomart query containing Ensembl gene ID, gene start and end position, chromosome, gene name (symbol), and strand
- `[OUT_STEM].locations.bed`: BED file containing all the locations of the genes selected for inclusion in the element set
- `[OUT_STEM].intersections`: Calculated intersections/bp between HARs for randomly generated controls and provided element set

### sigtest
Inputs  
- `[INPUT]`: Filepath to results file from HARE with simulation and element set intersection/bp values. Can also be a comma-separated list

Outputs  
- `[OUT_STEM].stats`: Tab-separated text file with the test parameters and resulting p-value
- `[OUT_STEM].png`: Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set

### Arguments and Options
You and also view this using the `-h` option after installation.

#### intersect
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--gwas`, `-g` | string | (Required) Filepath for GWAS summary statistics. |
| `--ref` | string | (Required) Filepath for BED-format human genome reference annotation. This file is used to simulate matched element sets. UCSC GRCh37 and GRCh38 gene annotation files are available under `reference_assets`. |
| `--eoi` | string | Filepath for BED-format elements of interest (EOIs). Required except when using `--anno_only` option (see below). The Richard, et al., 2020 HAR supplement files (original GRCh37 and a GRCh38 liftover) are available in hare.reference.assets.tar.gz . |
| `--out`, `-o` | string | Output filepath stem. If nothing is provided, will use the stem from the input GWAS file. |
| `--pval`, `-p` | float | P-value threshold for inclusion of variants in the annotation. Default is 1e-6. |
| `--gwas_p` | string | Column name for p-values in GWAS summary statistics. Default is "P". |
| `--maf`, '-m' | float | Minimum minor allele frequency (MAF). Default is 0.01. |
| `--gwas_maf` | string | Column name for minor allele frequency (MAF) in GWAS summary statistics. Default is "MAF". |
| `--source_neale` | - | Use Neale Lab GWAS summary statistics format. Will take priority over --gwas_p and --gwas_maf options. |
| `--source_bolt` | - | Use BOLT-LMM GWAS Summary Statistics format. Will take priority over --gwas_p and --gwas_maf options. |
| `--anno_only` | - | Only perform annotation, do not simulate element sets and perform intersections. EOI file is not required when this option is used. |
| `--keep_tmp` | - | Keep temporary files. Can be useful for troubleshooting. |
|`--ref_build`, `-r` | string | You may provide a human genome reference build which is used for the annotation. Options are either `37` (for GRCh37, hg19) or `38` (for GRCh38, hg38). Default value is `37`. |
|`--dist`, `-d` | int | Distance to transcript for which VEP assigns upstream and downstream consequences. Default is 5,000 bp. For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#output). |
|`--cache_dir` | string | VEP cache directory to use. Default is "$HOME/.vep/". For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#cacheopt) |
|`--cache_version` | string | VEP cache version to use. Default is 105 (version used during development).|
|`--biotypes` | string | Allowed biotypes for annotation. Options are \"protein_coding\", \"protein_all\", and \"all_features\". Default is \"protein_all.\" See [VEP's biotype documentation](https://uswest.ensembl.org/info/genome/genebuild/biotypes.html) for details.|
| `--draws`, `-n` | int | Number of simulations (draws) to use for background distribution. Default is n=1000. |

#### sigtest
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--input`, '-i' | string | Filepath for output of hare intersect (by default will have name [FILE].intersections). |
| `--out`, `-o` | string | Output filepath stem. If nothing is provided, will use "hare". |
| `--distribution` | str | Provide a distribution to use for hypothesis testing. By default will provide empirical p-values as well as those derived from Weibull distribution. Options are WEIBULL, NORMAL, BETA, and GAMMA (case insensitive). |
| `--skip_plot` | - | Do not generate plots. Default is OFF (will plot by default). |

## Reporting Bugs and Issues
To report bugs, please create an issue here or contact me via osmith@utexas.edu.

## Citation
If you use this pipeline in your work, please cite our article:

Kun E, Javan EM, Smith OS, et al., The genetic architecture and evolution of the human skeletal form.
Science 381,eadf8009(2023). DOI:10.1126/science.adf8009
