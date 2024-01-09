# HARE File Descriptions and Data Dictionaries
Olivia Smith  
osmith@utexas.edu  
GitHub: @ossmith

Below is comprehensive information about each of the files created using the HARE pipeline and data dictionaries for those which require it.

## hare.reference.assets.tar.gz
### HAR BED Files
`harsRichard2020.GRCh37.bed` and `harsRichard2020.GRCh38.bed`
BED file which contains the human accelerated regions (HARs) discovered/annotated in various different publications. This file is sourced from the supplement of [Richard, et al., 2020](https://doi.org/10.1016/j.cell.2020.02.057). The original HAR file was created with using human genome reference GRCh37. The GRCh38 HAR BED file was created using the UCSC Genome Browser liftover tool. Any file in BED format can be used to build the elements of interest set.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome |
| CHR_START | int | Starting position of the HAR (bp) |
| CHR_END  | int | End position of the HAR (bp) |
| PUB | string | Publication which identified/annotated the HAR |

### UCSC Genome Annotations
`UCSC.GRCh37.annotation.autosomes.bed.gz` and `UCSC.GRCh37.annotation.autosomes.bed.gz` are gzipped BED files with all of the gene annotations for human autosomes in GRCh37 and GRCh38. Each autosome was downloaded separately from [UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19) and then assembled into the single autosome file. Unzip this file for viewing or for use in the pipeline with the following command:
```
gunzip UCSC.GRCh37.autosomes.bed.gz
```
Data dictionary descriptions are taken from the UCSC Genome Browser at time of publication.

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

## [OUT_STEM].snps
A list of SNPs with genome-wide association to the phenotype in VCF file format. See [SAMtools documentation](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for more details on this file format.

## [OUT_STEM].annotation
Output from the Ensembl Variant Effect Predictor command line tool. Comments are preceded by `#`, including column headers. For data dictionary, see the ['default VEP output documentation'](https://uswest.ensembl.org/info/docs/tools/vep/vep_formats.html#defaultout).

## [OUT_STEM].biomart
Output from the BioMart location finding for the elements annotated by VEP (`[OUT_STEM].annotation` file). Headers are included in this file.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| ENSEMBL_ID | string | EnsemblID for the feature |
| START | int | First position of feature |
| END  | int | Last position of feature |
| CHR | int | Chromosome the feature is located on |
| GENE_NAME | string | Gene name associated with the feature in Ensembl |
| STRAND | integer | Strand which the feature is located on. `1` is for forward and `-1` is for reverse. |

## [OUT_STEM].locations.bed
BED file which contains only the locations of the elements annotated via VEP which will be intersected against the genomic elements of interest.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| CHR | int | Chromosome |
| START | int | Starting position of the feature (bp) |
| END  | int | End position of the feature (bp) |

## [OUT_STEM].intersections
This file contains the calculations of the intersections/bp for the simulation and phenotype-associated element sets.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| category | string | Category for the calculation which specifies whether the element set is either a `simulation` or the `test_set` (phenotype-associated). |
| int_per_bp | float | Intersections per base pair computed across the entire element set. |
| set_size  | int | Number of elements present in the element set. This number should be the same across all simulations associated with a given phenotype-associated element set. |

## [OUT_STEM].stats
sigtest results file which contains information about the run parameters, model fitting, and hypothesis testing (including the p-value) of the intersect results.

| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| FILENAME | string | Category for the calculation which specifies whether the element set is either "simulation" or "test_set" (phenotype-associated). |
| SET_SIZE | int | Number of elements present in the element set. This number applies both to the phenotype-associated and simulation element sets. |
| N_SIMULATIONS  | int | Number of simulations used to generate the background distribution. |
| SIM_IPB | float | The mean intersections/bp across all the simulations. |
| SET_IPB | float | The intersections/bp in the phenotype-associated element set. |
| P_EMPIRICAL | float | Empirical p-value calculated as the fraction of the simulations which had higher intersections/bp than the phenotype-associated element set. |
| WEIBULL_SHAPE | float | INFO |
| WEIBULL_SCALE | float | INFO |
| P_WEIBULL | float | P-value calculated from fit to weibull distribution (one tailed). |
| ADJUSTED_P | float | Empirical p-value adjusted for multiple hypothesis testing using Benjamini-Hochberg method. Note that this value should only be considered valid if all tests are provided as input in one run (use comma-separated list, see README for details). |

## [OUT_STEM].rnk
prerank results file which contains ranked list of genes and a score (either minimum or mean depending on score method used). This file does not contain headers.

| Column | Data Type | Description |
| ----------- | --------- | ----------- |
| 1 | string | HGNC symbol for the feature or gene |
| 2 | float | Associated score. Computed as the minimum or mean -log10(p) of positions associated with that feature. |
