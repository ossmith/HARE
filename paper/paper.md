---
title: 'HARE: A Python workflow for analyzing genomic feature enrichment in GWAS datasets'
tags:
  - Python
  - bioinformatics
  - evolution
  - gwas
  - enrichment analysis
authors:
  - name: Olivia S. Smith
    orcid: 0000-0001-5435-2982
    affiliation: "1"
  - name: Eucharist Kun
    orcid: 0000-0001-9919-173X
    affiliation: "1"
  - name: Vagheesh M. Narasimhan
    orcid: 0000-0001-8651-8844
    affiliation: "1, 2, 3"
affiliations:
  - name: Department of Integrative Biology, The University of Texas at Austin, United States of America
    index: 1
  - name: Department of Statistics and Data Science, The University of Texas at Austin, United States of America
    index: 2
  - name: Department of Population Health, Dell Medical School
    index: 3
date: 2 January 2024
bibliography: paper.bib
---

# Summary
Genomic datasets have grown by orders of magnitude in the last decade and thus require increasingly flexible and streamlined analysis workflows. One major line of research in genomics is gene enrichment analysis, which aims to quantitatively assess whether the genomic basis of a phenotype of interest has higher association with other genomic regions than would be expected by chance. HARE is a Python pipeline which annotates genome-wide significant positions, identifies overlap between these positions and a provided list of genomic features (here called 'elements') of interest, and determines enrichment likelihood through a resampling process of random length-matched regions. HARE is written in Python and is available for installation from source or via the Python Package Index (pip).

# Statement of need
HARE is a computational pipeline run via command line for analyzing enrichment of element sets, such as the set of human-accelerated regions (HARs), in genomic regions associated with a phenotype of interest. Several genetic enrichment approaches are already available in the literature, such as gene set enrichment analysis (GSEA) [@Subramanian2005; @Mootha2003], stratified linkage disequilibrium score regression (LDSC) [@Finucane2015], and FUMA GWAS [@Wantanabe2017]. However, these approaches are not translatable to studies examining overlaps in custom regions of the genome or when these regions are short, which is known to lead to bias in some approaches [@Finucane2015]. Our pipeline is flexible on either end of the enrichment analysis, both in terms of the regions it utilizes as well as the type of input genomic scan provided, allowing for a completely flexible model for genomic enrichment analysis.

Many tools also require time-consuming data reformatting which HARE minimizes by accepting a variety of summary statistic file formats. HARE dependencies can be installed with conda and pip, allowing for an immediate run of the complete analysis pipeline. HARE has already been used to analyze over 200 traits in published work as well as ongoing studies presented at conferences [@Kun2023a; @Kun2023b; @Xu2023], exemplifying its usefulness in current and future research.

# Dependencies
HARE requires Ensembl's Variant Effect Predictor [@McLaren2016] with the appropriate cache ([VEP Cache Documentation](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)) as well as BEDTools [@Quinlan2010]. It also requires the use of `wget` for querying Ensembl BioMart [@Martin2023]. It is recommended that users install dependencies through a combination of conda using the `environment.yml` file included in the repository and pip. Complete instructions for dependency and tool installation can be found in the repository README.

# Workflow
A visual depiction of the HARE workflow is shown in Figure 1.

![HARE visual workflow. **(A)** Summary statistics from a genome-wide association study (GWAS) are filtered to genome-wide significant single nucleotide polymorphisms (SNPs) based on p-value threshold. These are the phenotype-associated (PA) SNPs. **(B)** For each PA SNP, the nearest genetic feature within a given bp distance and its location are annotated. This is called the element set. **(C)** Simulated element sets built from random regions matched on length to the PA features are generated. **(D)** Intersections per bp between each element set, both PA and simulated, and the elements of interest (e.g. human accelerated regions) and an empirical p-value are computed.](Fig1_HAREWorkflow.png "HARE Workflow")

## Annotation
The first portion of the workflow involves identifying and annotating the phenotype-associated (PA) features to be tested for enrichment (Figure 1A-B). First, genome-wide significant single nucleotide polymorphisms (SNPs) from the GWAS summary statistics are identified using a p-value threshold (e.g. $p < 1\times10^{-8}$). Using Ensembl's Variant Effect Predictor (VEP), the nearest gene within an upstream/downstream buffer distance (bp) is annotated. Users can specify what annotation biotypes are allowed -- protein coding, all protein biotypes, or regulatory features. We then use Ensembl's BioMart to locate the chromosome, start, and stop positions of all of these features. These features compose the phenotype-associated (PA) element set for enrichment analysis.

## Simulation-based enrichment analysis
Once the element set is defined, we determine parameters for creation of a matched set of genetic regions in the genome (Figure 1C). We split the element set into equally-sized bins and compute the average length within each bin. Using these parameters, we simulate a large number of element sets with BEDTools random [@Quinlan2010] which then comprise a background distribution for testing enrichment. We recommend that users use a minimum of 1,000 simulations and element sets of at least 40 features for analysis. With BEDTools intersect, we calculate the intersections per base pair of sequence between the phenotype-associated features and the comparison set (the elements of interest, e.g. HARs). To determine enrichment, we compute an empirically-computed p-value comparing the observed intersections/bp with the distribution obtained from the random re-sampling (one tailed, Figure 1D).

## Prerank
For additional enrichment analysis, tools such as GSEA [@Subramanian2005; @Mootha2003] and WebGestalt [@Zhang2005; @Wang2013; @Wang2017; @Liao2019] require ranked lists of genes. HARE is able to create this ranked list using any genome-wide summary statistic (typically results of GWAS or selection scans) containing genomic positions and associated p-values. It annotates genes within an upstream/downstream buffer distance (filtering on top hits or based on p-value threshold is allowed) and builds a score for each gene. This score is the average $-log_{10}$(p-value) for all positions associated with the gene. This gene:score ranked list is provided and can be used with these additional enrichment analysis tools.

# Example case
In @Kun2023a, we used HARE to analyze enrichment of human accelerated regions (HARs) within gene sets associated with skeletal phenotypes. We performed automated image processing of DXA x-rays from over 30,000 individuals from the UK Biobank to generate measurements of bone lengths. We then performed genome wide association studies (GWAS) for 23 image-derived phenotypes such as hip width to shoulder width, arm to leg, and tibia to femur ratios. Using HARE, we identified genome-wide significant SNPs associated with each phenotype as well as a number of other phenotypes with publicly-available GWAS summary statistics and created phenotype-associated gene sets. Through comparing these gene sets to human accelerated regions, we found enrichment in 8 phenotypes (Figure 2).

![HARE example cases showing p-values of enrichment for overlap between skeletal, dermatological, endocrine, neurological, cancer, metabolic, autoimmune, and gastrointestinal phenotypes and human accelerated regions (HARs) as compared to randomly sampled gene sets of comparable length distribution. Traits with FDR-corrected p-values of less than 0.05 are shown in orange and traits above the threshold are shown in blue.](Fig2_HAREExampleCases.png "HARE Example Cases")

We have also tested HARE for use in non-human species based on GWAS datasets from *Bos taurus* [@Zhuang2020] (vertebrate), *Arabidopsis thaliana* [@Atwell2010] (plant), *Drosophila melanogaster* and *Apis mellifera* [@Avalos2020] (metazoa), and *Saccharomyces cerevisiae* [@Sardi2018] (fungus).

# Acknowledgements
The authors would like to thank Liaoyi Xu for providing additional user testing and feedback. O.S.S. was supported by NSF Graduate Research Fellowship (DGE 2137420). GPU and compute resources were supported by a Director’s Discretionary Award from the Texas Advanced Computing Center.

# References
paper.bib
