# HARE Genetic Feature Enrichment Analysis Pipeline
Author: Olivia Smith, osmith@utexas.edu, GitHub [@ossmith](https://github.com/ossmith)  
Publication Corresponding Author: Vagheesh Narasimhan, vagheesh@utexas.edu, GitHub [@vagheesh](https://github.com/vagheesh)

This script identifies the genes and genomic features associated with autosomal genome-wide significant positions based on GWAS summary statistics for a given phenotype. It then tests for elevated intersections between the phenotype-associated features and the elements of interest (e.g. HARs) based on a simulated background distribution.

## Repository structure
```
HARE
|
│-- README.md
|-- dataDictionary.md
|-- environment.yml
|-- pyproject.toml
|-- MANIFEST.in
|-- hare.reference.assets.tar.gz
|
|---src
│   │-- __init__.py
│   │-- hare.py
│   │-- hareclasses.py
│   │-- intersect.py
│   │-- sigtest.py
│   │-- prerank.py
|
|---example
│   │-- exampleREADME.md
│   │-- example_run.sh
|   |-- example_{filename.extension}
│
|---tests
│   │-- __init__.py
│   │-- intersect_test.py
│   │-- sigtest_test.py
│   │-- prerank_test.py
│   │-- species_test.py
│   |-- input
|       |-- ...
```

## Installation
### Dependencies
The following are required dependencies for HARE. Hyperlinks will direct you to the installation documentation for these dependencies. You may either install these independently or through creation of a conda environment using the provided `environment.yml` file (see below).

- Linux OS (<sup>*</sup>see note below regarding MacOS/Unix compatibility)
- Python => 3.8, pip, and the following packages: argparse, matplotlib, numpy, pandas, scipy>1.8
- [VEP command line tool](https://ensembl.org/info/docs/tools/vep/script/vep_download.html) with [reference genome cache](https://ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) (for non-human species, find the caches at their respective Ensembl sister sites)
- perl>5.32 and [perl Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)
- [HTSlib](http://www.htslib.org/download/)
- [BEDtools v2.30.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- wget

<sup>*</sup>Note: HARE has been tested and is functional in Unix/MacOS environments where dependencies can be installed, but due to some MacOS incompatibilities with Ensembl VEP, it may not be possible to install the necessary dependencies using the `environment.yml` file. Users can see detailed instructions for installation from Ensembl [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#macos) and follow links above for the other dependencies.

### Install with conda and pip
Start by cloning this repo:
```
git clone https://github.com/ossmith/HARE.git
```

You can install the dependencies as a conda environment using the provided `environment.yml` file. By default it will create an environment named 'hare-env'.  
```
cd HARE
conda env create -f environment.yml
conda activate hare-env
```

You may then need to update the perl Compress::Raw::Zlib module:
```
cpan Compress::Raw::Zlib
```

And then build HARE from the included setup files with:
```
pip3 install .
```

**You will then need to download a VEP genome cache if you have not already done so.** Find the right cache and install using the instructions [here](https://ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).  

Example (*H. sapiens* GRCh37 VEP cache v105):
```
cd $HOME/.vep
wget https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_vep_105_GRCh37.tar.gz
```

Ensembl recommends that you use the same cache version as your VEP install version (hare-env includes v105). Note that the first time you run the tool after downloading the cache, it will need to do some configuration that will slow down the pipeline. This will only take place the first time you use a given cache.

Included reference assets can be extracted with:
```
tar xvf hare.reference.assets.tar.gz
```

### Testing installation
You can run the unit tests to confirm functionality (requires path to installed cache and version if different from defaults, defaults are '$HOME/.vep/' and '105') using:
```
cd tests/
python intersect_test.py --cache_dir {VEP_CACHE_PATH} --cache_ver {VEP_VERSION}
python sigtest_test.py
python prerank_test.py
```

Please note that while a test suite for functionality on non-human species is available (`species_test.py`), it requires all 5 non-human example species caches. We recommend instead that users interested in this function test only the species of interest with the datasets found under `tests/input/{ref/gwas/eoi}_{species}.txt` if necessary.

### Running an example
There is an example fileset in the `HARE/example` directory, including (i) example inputs, (ii) a bash script (`example_run.sh`) to run each command on the example files, and (iii) example output files. Note that this script still requires users to provide the path to their genome reference (the corresponding GRCh37 reference for this example is located in `hare.reference.assets.tar.gz`) and VEP cache. Refer to the included `exampleREADME.md` for additional file details.

## Commands
### Running annotation, simulation, and intersection
The primary analysis in this repository is performed with the Python `intersect` function. It performs annotation, simulation, and intersection for the provided GWAS and elements of interest (for intersecting against). The most basic command is:
```
hare intersect --gwas [GWAS] --eoi [EOI] --ref [REF] --out [OUT_STEM] ... [OPTIONS]
```  
Various options are explained below.

### Running significance testing
The significance testing is done using `sigtest`. This function is separated from the main analysis because it allows for multiple HARE output files from different GWASs to be analyzed together to create a condensed TSV. The command is:
```
hare sigtest --i [INPUT] --o [OUT_STEM] ... [OPTIONS]
```
where `--input` is the `[FILE].intersections` file (or comma-separated list of these files) created using `intersect`.

### Running gene preranking
Other enrichment analysis tools such as WebGestalt and GSEA make use of ranked gene lists (HGNC symbol:score pairings). HARE can annotate and rank genes from any dataset with genomic position and p-value data (e.g. selection scan, GWAS). It can also be used for peak calling by identifying genes with multiple positions of genome-wide significance (by p-value).
```
hare prerank --i [INPUT] --o [OUT_STEM] ... [OPTIONS]
```

## File Inputs and Outputs
### intersect
Inputs  
- `[GWAS]`: Tab-separated GWAS summary statistics file which contains information on the chromosome, position, MAF, and p-value for each SNP
- `[EOI]`: BED file with elements of interest to intersect against (EOIs). The Richard, et al., 2020 HARs BED file can be found in `hare.reference.assets.tar.gz`
- `[REF]`: Reference genome file containing chromosomes and associated lengths (see more details of what is required and how to make one [here](https://bedtools.readthedocs.io/en/latest/content/general-usage.html?highlight=genome%20file#genome-file-format)). Reference genome files for gene annotation of GRCh37 and GRCh38 can be found in `hare.reference.assets.tar.gz`

Outputs  
- `[OUT_STEM].snps`: List of variants which pass QC conditions (biallelic, MAF threshold, p-value threshold) and will be used for annotation and analysis
- `[OUT_STEM].annotation`: Raw VEP annotation results
- `[OUT_STEM].annotation_summary.html`: VEP run log produced by command line run
- `[OUT_STEM].biomart`: Output from Biomart query containing Ensembl gene ID, gene start and end position, chromosome, gene name (symbol), and strand
- `[OUT_STEM].locations.bed`: BED file containing all the locations of the genes selected for inclusion in the element set
- `[OUT_STEM].intersections`: Calculated intersections/bp between EOIs (e.g. HARs) and simulated or phenotype-associated element sets

### sigtest
Inputs  
- `[INPUT].intersections`: Filepath to results file from `intersect` (naming convention `[OUT_STEM].intersections`) with simulation and element set intersection/bp values. Can also be a comma-separated list of such filenames

Outputs  
- `[OUT_STEM].stats`: Tab-separated text file with the test statistics and resulting p-value
- `[OUT_STEM].png`: Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set

### prerank
Inputs
- `[INPUT]`: File containing p-values associated with each position. Chromosome, position, and p-value columns are required.
- `[INCLUDE].bed` (optional): BED format file containing a list of positions to filter to for analysis. Only listed positions will be included.
- `[EXCLUDE].bed` (optional): BED format file containing a list of positions to remove from analysis.

Outputs
- `[OUTPUT].rnk`: File containing HGNC symbol and score for each annotated gene. Use `--score_method` to choose whether min or mean p-value from associated positions is used to score.
- `[OUTPUT].dmp` (if `--call_peaks` used): File containing peaks for which a sufficient number of positions within a given distance meet genome-wide significance threshold. Example: peaks with at least 3 positions within 1,000 bp of a gene have a p-value below 1e-06.

### Arguments and Options
You can also view this using the `hare {FUNCTION} -h` option after installation.

#### intersect
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--gwas`, `-g` | string | (Required) Filepath for GWAS summary statistics. |
| `--ref` | string | (Required) Filepath for genome reference annotation (either a BED or fai file works). This file is used to simulate matched element sets. Find more details about what is required [here](https://bedtools.readthedocs.io/en/latest/content/general-usage.html?highlight=genome%20file#genome-file-format) GRCh37 and GRCh38 genome files are available in `hare.reference.assets.tar.gz`. |
| `--eoi` | string | Filepath for BED-format elements of interest (EOIs). Required except when using `--anno_only` option (see below). The Richard, et al., 2020 HAR supplement files (original GRCh37 and a GRCh38 liftover) are available in `hare.reference.assets.tar.gz`. |
| `--out`, `-o` | string | Output filepath stem. If nothing is provided, will use the stem from the input GWAS file. |
| `--pval`, `-p` | float | P-value threshold for inclusion of variants in the annotation. Default is 1e-6. |
| `--gwas_p` | string | Column name for p-values in GWAS summary statistics. Default is "P". |
| `--use_z` | - | Use Z score (in "Z" column) to calculate p-values for GWAS summary statistics. By default is OFF. |
| `--maf`, `-m` | float | Minimum minor allele frequency (MAF). Default is 0.01. |
| `--gwas_maf` | string | Column name for minor allele frequency (MAF) in GWAS summary statistics. Default is "MAF". |
| `--gwas_ref` | string | Column name for reference allele in GWAS summary statistics. Default is "REF". |
| `--gwas_alt` | string | Column name for alternate allele in GWAS summary statistics. Default is "ALT". |
| `--snp_map` | string | If SNPs provided as IDs instead of genomic locations (CHR, POS), provide a BED file which maps IDs to locations. |
| `--source_neale` | - | Use Neale Lab GWAS summary statistics format. Will take priority over `--gwas_{COLUMN_NAME}` options. |
| `--source_bolt` | - | Use BOLT-LMM GWAS Summary Statistics format. Will take priority over `--gwas_{COLUMN_NAME}` options. **Uses `P_BOLT_LMM_INF` for p-values.** |
| `--anno_only` | - | Only perform annotation, do not simulate element sets and perform intersections. EOI file is not required when this option is used. |
|`--ref_build`, `-r` | string | Genome reference build (only for human data). Options are either 37 (for GRCh37, hg19) or 38 (for GRCh38, hg38). Default is 37. |
|`--dist`, `-d` | int | Distance to transcript for which VEP assigns upstream and downstream consequences. Default is 5,000 bp. For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#output). |
|`--cache_dir` | string | VEP cache directory to use. Default is "$HOME/.vep/". For more details, see [VEP CLI documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_options.html#cacheopt) |
|`--cache_version` | string | VEP cache version to use. Default is 105 (version used during development). |
|`--species` | string | Default is "homo_sapiens". Latin name of the species for your data e.g. \"arabidopsis_thaliana\". If not human, pair with a sister flag (below, e.g. --plant, --metazoa). This is **required** in order to find the correct VEP cache and Ensembl site. |
|`--vertebrate` | string | Indicate data is from non-human, vertebrate species. Requires a corresponding species name using the --species flag. |
|`--plant` | string | Indicate data is from a plant species. Requires a corresponding species name using the --species flag. |
|`--metazoa` | string | Indicate data is from a metazoic species. Requires a corresponding species name using the --species flag. |
|`--fungi` | string | Indicate data is from a fungal species. Requires a corresponding species name using the --species flag. |
|`--biotypes` | string | Allowed biotypes for annotation. Options are \"protein_coding\", \"protein_all\", and \"all_features\". Default is \"protein_all.\" See [VEP's biotype documentation](https://uswest.ensembl.org/info/genome/genebuild/biotypes.html) for details. |
| `--draws`, `-n` | int | Number of simulations (draws) to use for background distribution. Default is n=1000. |
| `--keep_tmp` | - | Keep all temporary files created during the run. |

#### sigtest
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--input`, `-i` | string | (Required) Filepath for output of `intersect` (default output has file extension `.intersections`). |
| `--out`, `-o` | string | Output filepath stem. If nothing is provided, will use "hare". |
| `--skip_plot` | - | Do not generate plots. Default is OFF (will plot by default). |

#### prerank
| Option | Data Type | Description |
| ----------- | --------- | ----------- |
| `--input`, `-i` | string | (Required) Filepath of file with p-values. |
| `--output`, `-o` | string | Output filepath stem. If nothing is provided, will use input stem. |
| `--ref_build`, `-b` | string | GRCh reference build. Options are either 37 (hg19) or 38 (hg38). Default is 37. |
| `--topN`, `-n` | int | Number of lowest p-value positions to use. No default (all positions will be scored and ranked). |
| `--pval`, `-p` | float | P-value threshold for annotated genes. Default is 1 (all positions will be scored and ranked). |
| `--pval_col`, `-vc` | string | Column with p-values. Default is 'P'. |
| `--chr`, `-c` | int | Filter to a single given chromosome (multi-chromosome or chromosome:region syntax is not supported). Will analyze all positions if not provided. Recommended for large (whole genome) datasets. |
| `--chr_col`, `-cc` | string | Column with chromosome. Default is 'CHR'. |
| `--pos_col`, `-pc` | string | Column with position. Default is 'POS'. |
| `--incl`, `-f` | string | Filepath of BED file with positions to include in analysis. |
| `--excl`, `-x` | string | Filepath of BED file with positions to exclude from analysis. |
| `--buffer` | int | Number of bp to look upstream and downstream for annotations. Default is 0 (will only annotate genes if it is within the bounds). |
| `--biotypes`, `-t` | str | Allowed gene types when annotating. Options are 'all' or 'coding'. Default 'all'. See HGNC documentation for gene type details. |
| `--score_method`, `-m` | str | Method for how to compute score. Options are 'mean' or 'min'. Mean will use the average p-value of all positions annotated to the gene and min will select the minimum p-value of those positions. Default is minimum. |
| `--call_peaks` | - | Call peaks. See --dmp_{OPTION} options for settings. Default is OFF (does not call peaks). |
| `--dmpN` | int | Number of positions annotating to a given gene required within provided distance for peak calling. Default is 3. Set to 1 to retrieve all. |
| `--dmpP` | float | P-value threshold for peak (differential methylation) calling. Default is 1e-08. |
| `--dmpD` | int | Distance (bp) allowed for which dmpN and dmpP conditions must be met for peak calling. If none provided, will use --buffer distance. |

## Reporting Bugs and Contribution Guidelines
Please submit issues, bug reports, and feature requests using GitHub Issues.  

If you would like to contribute to the code, please fork this repository and then submit a pull request (PR) with your changes. We require that any changes or additions to functionality also include a corresponding automated test.

## Citation
If you use this pipeline in your work, please cite these articles:

Smith OS, Kun E, and Narasimhan VM, (2024). HARE: A Python workflow for analyzing genomic feature enrichment in GWAS datasets. Journal of Open Source Software, 9(97), 6359, https://doi.org/10.21105/joss.06359

Kun E, Javan EM, Smith OS, et al., (2023). The genetic architecture and evolution of the human skeletal form.
Science 381, eadf8009, https://doi.org/10.1126/science.adf8009
