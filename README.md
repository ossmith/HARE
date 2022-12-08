# Human Accelerated Region Enrichment (HARE) Evolutionary Analysis Pipeline
Author: Olivia Smith, osmith@utexas.edu, GitHub @ossmith  
Publication Corresponding Author: Vagheesh Narasimhan, vagheesh@austin.utexas.edu, GitHub @vagheesh

This script identifies the genes and genomic features associated with autosomal genome-wide significant SNPs based on GWAS summary statistics for a given phenotype. It then tests for elevated intersections between the genomic features and human accelerated regions (HARs) based on a simulated background distribution. This analysis is based upon a modification of the method presented in [Richard, et al., 2020](https://doi.org/10.1016/j.cell.2020.02.057).

## Repository Structure
```
HARE
│   README.md    
|   dataDictionary.md
|   HARE.py
|   HARE_Results.R
|   environment.yml
│
└───reference_assets
│   │   harsRichard2020.GRCh37.bed
│   │   harsRichard2020.GRCh38.bed
│   │   UCSC.GRCh37.autosomes.annotation.bed.gz
│   │   UCSC.GRCh38.autosomes.annotation.bed.gz
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
The following are required dependencies for HARE. Hyperlinks will direct you to the installation documentation for these dependencies.
- Python => 3.0
- [VEP command line tool](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html) with [human reference genome cache](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)
- [perl Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)
- [HTSlib](http://www.htslib.org/download/)
- [BEDtools v2.30.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- wget (Note that wget is installed by default on Linux but not Mac (Unix) operating systems, so you may need to install it separately, e.g. `brew install wget`.)
- [R => 4.1](https://www.r-project.org/) (only required for `HARE_Results.R`)

You can also install the dependencies as a conda environment using the provided `environment.yml` file. Once in the HARE directory, create and activate the conda environment using:
```
conda env create -f environment.yml
conda activate HARE
```

You may then need to update the perl Compress::Raw::Zlib module:
```
cpan Compress::Raw::Zlib
```

### Installation from source
Clone this repository:
```
git clone https://github.com/ossmith/HARE.git
```

## Workflow
### Inputs
- `[GWAS]`: Tab-separated GWAS summary statistics file which contains information on the chromosome, position, MAF, and p-value for each SNP
- `[HAR]`: BED file with human accelerated regions (HARs). The Richard, et al., 2020 HARs BED file can be found in the `reference_assets` directory.
- `[REF]`: BED file with gene annotation of human genome reference (e.g. UCSC's hg19 .bed). UCSC genome annotations for GRCh37 and GRCh38 can be found in the `reference_assets` directory.

### Outputs
- `[OUT_STEM].snps`: List of variants which pass QC conditions (biallelic, MAF threshold, p-value threshold) and will be used for annotation and analysis
- `[OUT_STEM].annotation`: Raw VEP annotation results
- `[OUT_STEM].annotation_summary.html`: VEP run log produced by command line run
- `[OUT_STEM].biomart`: Output from Biomart query containing Ensembl gene ID, gene start and end position, chromosome, gene name (symbol), and strand
- `[OUT_STEM].locations.bed`: BED file containing all the locations of the genes selected for inclusion in the element set
- `[OUT_STEM].intersections`: Calculated intersections/bp between HARs for randomly generated controls and provided element set
- `[OUT_STEM].results.tsv`: (HARE_Results.R) Tab-separated file with the test parameters and resulting p-value
- `[OUT_STEM].results.png`: (HARE_Results.R) Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set

## Commands
### Running HARE.py
The primary analysis in this repository is performed with the Python `HARE.py` script. The most basic command for this script is
```
python HARE.py --gwas [GWAS] --HAR [HAR] --ref [REF] --out [OUT_STEM] ... [OPTIONS]
```  
The various options are explained below.

### Running HARE_Results.R
The hypothesis testing is performed by the `HARE_Results.R` script. This file is separated from the main script because it allows for multiple HARE output files from different GWASs to be analyzed together to create a condensed TSV. The command is always:
```
R-script HARE_Results.R --i [INPUT] --o [OUT_STEM]
```
where `--input` is the `[FILE].intersections` file (or comma-separated list of these files) created by `HARE.py`.

### Options
#### Basic Arguments
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--gwas`, `-g` | string | (Required) Filepath for GWAS summary statistics. |
| `--ref` | string | (Required) Filepath for BED-format human genome reference annotation. This file is used to simulate matched element sets. UCSC GRCh37 and GRCh38 gene annotation files are available under `reference_assets`. |
| `--HAR` | string | (Required) Filepath for BED-format human accelerated regions (HARs). The Richard, et al., 2020 HAR supplement files (GRCh37 and a GRCh38 liftover) are available under `reference_assets`. |
| `--out`, `-o` | string | Output filepath stem. If nothing is provided, will use the stem from the input GWAS file. |

#### GWAS Options
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--pval`, `-p` | float | P-value threshold for inclusion of variants in the annotation. Default is 1e-6. |
| `--gwas_p` | string | Column name for p-values in GWAS summary statistics. Default is "P". |
| `--maf` | float | Minimum minor allele frequency (MAF). Default is 0.01. |
| `--gwas_maf` | string | Column name for minor allele frequency (MAF) in GWAS summary statistics. Default is "MAF". |
| `--source_neale` | string | Use Neale Lab GWAS summary statistics format. Will take priority over --gwas_p and --gwas_maf options. Accepts T, True, F, and False (case insensitive). Default False. |
| `--source_bolt` | string | Use Javan GWAS Summary Statistics format. Will take priority over --gwas_p and --gwas_maf options. Accepts T, True, F, and False (case insensitive). Default False. |

#### Annotation Options (Ensembl VEP and BioMart)
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
|`--ref_build`, `-r` | string | You may provide a human genome reference build which is used for the annotation. Options are either `37` (for GRCh37, hg19) or `38` (for GRCh38, hg38). Default value is `37`. |
|`--dist`, `-d` | int | Distance to transcript for which VEP assigns upstream and downstream consequences. Default is 5,000 bp. For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#output). |
|`--cache_dir` | string | VEP cache directory to use. Default is "$HOME/.vep/". For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#cacheopt) |
|`--cache_version` | string | VEP cache version to use. Default is 105 (version used during development).|
|`--biotypes` | string | Allowed biotypes for annotation. Options are \"protein_coding\", \"protein_all\", and \"all_features\". Default is \"protein_all.\" See [VEP's biotype documentation](https://uswest.ensembl.org/info/genome/genebuild/biotypes.html) for details.|

#### Simulation Options
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--draws`, `-n` | int | Number of simulations (draws) to use for background distribution. Default is n=1,000. |

## Reporting Bugs and Issues
To report bugs, please create an issue or contact me via osmith@utexas.edu.

## Citation
If you use this pipeline in your work, please cite our preprint.
