#!/usr/bin/env Rscript
# ptest.R
# Olivia Smith
# osmith@utexas.edu
# GitHub: @ossmith

list.of.packages <- c("optparse", "fitdistrplus", "goft", "actuar")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){ stop(paste0("FileNotFoundError: some required R packages are missing. Cannot find the following packages: "), new.packages) }

suppressMessages(library(optparse,quietly=T))
suppressMessages(library(fitdistrplus,quietly=T))
suppressMessages(library(goft,quietly=T))
suppressMessages(library(actuar,quietly=T))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Filepath to results file from HARE.py ([OUTPUT].intersections) with simulation and element set intersection/bp values.", metavar="character"),
  make_option(c("-d", "--dist", type="character", default="weibull", help="Distribution to use for hypothesis testing [default= %default]"))
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
inFile = toString(opt[1])
distrib = toString(opt[2])

hypothesisTest <- function(filename, dist){
  # File intake
  raw_df <- read.table(filename, sep="\t",header=T)

  # Separate test set results from simulations
  test_set <- subset(raw_df, raw_df$category=="test_set")
  simulation_set <- subset(raw_df, raw_df$category=="simulation")
  nsim = nrow(simulation_set)
  simulation_set$int_per_bp[simulation_set$int_per_bp == 0] <- 1e-80
  # set_size <- test_set[1,3]

  if(length(unique(simulation_set$int_per_bp))==1){
    cat(stop("Error: all simulated intersection/bp are equivalent.
    This may be a result of an insufficient element set (we recommend set_size<=40) or insufficient simulations (we recommend n<=1000).
    If this error persists, please open an issue on GitHub."))
  }

  if (dist == "normal"){
    dist = "norm"
  }

  # Fit simulations to distribution and perform hypothesis testing on test set results
  fit_mod <- fitdist(simulation_set$int_per_bp, dist, lower=c(0,0))
  stat <- gofstat(fit_mod, fitnames=c(dist))

  if(dist == "weibull"){
    pval <- pweibull(test_set$int_per_bp, shape=fit_mod$estimate["shape"], scale=fit_mod$estimate["scale"], log.p = FALSE, lower.tail=FALSE)
  }
  else if(dist == "beta"){
    pval <- pbeta(test_set$int_per_bp, shape1=fit_mod$estimate["shape1"], shape2=fit_mod$estimate["shape2"], log.p = FALSE, lower.tail=FALSE)
  }
  else if(dist == "gamma"){
    pval <- pgamma(test_set$int_per_bp, shape=fit_mod$estimate["shape"], rate=fit_mod$estimate["rate"], log.p = FALSE, lower.tail=FALSE)
  }
  else if(dist == "norm"){
    pval <- pnorm(test_set$int_per_bp, mean=fit_mod$estimate["mean"], sd=fit_mod$estimate["sd"], log.p = FALSE, lower.tail=FALSE)
  }

  # Determine empirical p-value and return results
  fract_more <- length(simulation_set$int_per_bp[simulation_set$int_per_bp > test_set$int_per_bp])/length(simulation_set$int_per_bp)
  row <- c(nsim, dist, stat$ks[dist], mean(simulation_set$int_per_bp), test_set$int_per_bp, pval, fract_more)
  return(row)
}

# Print result to standard out (will be read in by Python script)
cat(noquote(hypothesisTest(inFile, distrib)))
