###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: check blank of extraction
#PROJECT: marsh_metabarcoding
###.............................................................................
library(phyloseq)
library(dplyr)
library(decontam)

#read phyloseq objects
ps <- readRDS("data/intermediate/ps.rds")

#load function to plot prevalence~abundace per phylum
source("code/functions/phyla_abundance_phyloseq.r")

#create phyloseq object with the negative control
ps_neg_prevalence <-
  phyloseq::prune_samples("BPCR", ps) %>%
  plot_abundances()
report_neg <-
 ps_neg_prevalence[[1]] %>%
  dplyr::filter(TotalAbundance > 0)
write.csv(report_neg, "output/report_neg.csv")
