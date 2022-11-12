#!/usr/bin/env Rscript
# HARE_Results.R
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

################################################################################

# This script performs model fitting and hypothesis testing on results generated
# using the HARE.py script. It then generates a PNG figure showing the value of 
# the intersections/bp of the phenotype-associated element set against the 
# background of simulated element sets.

# Dependencies: R packages optparse, fitdistpllus, goft, actuar; R => 4.1.0

# Inputs: [HARE_Results].intersections: Filename for HARE output (can also be comma-separated list) with calculated intersections/bp between HARs for randomly generated controls and provided element set

# Outputs: [OUT_STEM].tsv: Tab-separated file with the test parameters and resulting p-value
#          [OUT_STEM].png: Figure showing the distribution of intersections/bp of the simulations against the intersections/bp of the phenotype-associated element set


# Example command: Rscript HARE_Results.R --input [INPUT] --output [OUT_STEM]

################################################################################

script_version <- "1.0.0"

cat(noquote(paste0("\n\U0001F955-----Results processing started at ",Sys.time(),".\n")))

#### Load libraries and process inputs ####
cat(noquote("\nLoading libraries.\n"))
suppressMessages(library(optparse,quietly=T))
suppressMessages(library(tidyverse,quietly=T))
suppressMessages(library(fitdistrplus,quietly=T))
suppressMessages(library(goft,quietly=T))
suppressMessages(library(actuar,quietly=T))
cat(noquote("\nLibraries loaded."))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Filepath to results file from HARE.py ([OUTPUT].intersections) with simulation and element set intersection/bp values.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="HARE_Results.txt", 
              help="Output filepath [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
inFile = toString(opt[1])
inList = str_split(inFile,",")[[1]]
outFile = toString(opt[2])
cat(noquote(paste0("\nInput files: ",inFile)))

parseFile <- function(filename){
  # Load file
  cat(noquote("\nParsing files..."))
  raw_df <- read.table(filename, sep="\t",header=T)
  cat(noquote("OK"))
  return(raw_df)
}

hypothesisTest <- function(filename, raw_df){
  # Separate test set results from simulations
  test_set <- subset(raw_df, raw_df$category=="test_set")
  simulation_set <- subset(raw_df, raw_df$category=="simulation")
  nsim = nrow(raw_df)-1
  simulation_set$int_per_bp[simulation_set$int_per_bp == 0] <- 1e-80
  set_size <- test_set[1,3]
  
  # Fit simulations to a Weibull distribution and perform hypothesis testing on test set results
  cat(noquote("\nPerforming fitting and p testing..."))
  fit_wb <- fitdist(simulation_set$int_per_bp, "weibull",lower=c(0,0))
  stat <- gofstat(fit_wb, fitnames=c("weibull"))
  pval <- pweibull(test_set$int_per_bp, shape=fit_wb$estimate["shape"], scale=fit_wb$estimate["scale"], log.p = FALSE, lower.tail=FALSE)
  fract_more <- length(simulation_set$int_per_bp[simulation_set$int_per_bp > test_set$int_per_bp])/length(simulation_set$int_per_bp)
  fract_less <- length(simulation_set$int_per_bp[simulation_set$int_per_bp < test_set$int_per_bp])/length(simulation_set$int_per_bp)
  cat(noquote("OK"))
  row <- data.frame(file=filename, set_size=set_size, n_simulations=nsim, distribution="weibull", KS_stat=stat$ks["weibull"], mean_int_per_bp=mean(simulation_set$int_per_bp), set_int_per_bp=test_set$int_per_bp, pvalue=pval, frac_less=fract_less, frac_more=fract_more)
  return(list("row"=row,"p"=pval))
}

plotResults <- function(filename, raw_df, pval){
  # Create and save figure
  test_set <- subset(raw_df, raw_df$category=="test_set")
  fileExtless <- sapply(strsplit(filename,".intersections"), `[`, 1)
  out_png = paste0(fileExtless,".results.png")
  pval_txt = paste0("p-value: ",pval)
  binsize = (max(raw_df$int_per_bp)-min(raw_df$int_per_bp))/20
  harHist <- raw_df %>% subset(raw_df$category=="simulation") %>% ggplot() +
    aes(x=int_per_bp) +
    geom_histogram(binwidth=binsize, color="white", fill="#1F77B4") +
    ylab("Count") + xlab("Intersections/BP") +
    geom_vline(xintercept = test_set$int_per_bp, linetype="dotted", color=rgb(44,160,44,maxColorValue = 255), size=1.5) +
    annotate("text", Inf, -Inf, label=pval_txt, hjust = 1, vjust = -2) + theme_bw()
  ggsave(out_png, harHist, height=10, width=14)
}

#### Main ####
HEADER <- c("file", "set_size", "n_simulations", "distribution", "KS_stat", "mean_int_per_bp", "set_int_per_bp", "pvalue", "frac_less", "frac_more")
RESULTS.DF <- data.frame(matrix(ncol = length(HEADER), nrow=0))
colnames(RESULTS.DF) <- HEADER
for (i in inList){
  cat(noquote(paste0("\nAnalyzing ",i)))
  PARSED.DF <-  parseFile(i)
  TEST.OUT <- hypothesisTest(i, PARSED.DF)
  plotResults(i, PARSED.DF, TEST.OUT$p)
  RESULTS.DF <- RESULTS.DF %>% 
    rbind(TEST.OUT$row)
}

#### Write combined results table ####
cat(noquote(paste0("\nWriting results file to ",outFile,".results.tsv...")))
write.table(RESULTS.DF, paste0(outFile,".results.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
cat(noquote("OK"))

cat(noquote(paste0("\n\nResults processing completed at ",Sys.time(),".-----\U0001F407\n\n")))

