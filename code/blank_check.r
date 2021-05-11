###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: check blank of extraction
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(dplyr)
library(decontam)

#read phyloseq objects
ps_marsh <- readRDS("data/intermediate/ps_marsh.rds")
ps_greenh <- readRDS("data/intermediate/ps_greenh.rds")

#load function to plot prevalence~abundace per phylum
source("code/functions/phyla_abundance_phyloseq.r")

#create phyloseq object with the negative control
ps_neg_prevalence_marsh <-
  phyloseq::prune_samples("BPCR", ps_marsh) %>%
  plot_abundances()
report_neg_marsh <-
 ps_neg_prevalence_marsh[[1]] %>%
  dplyr::filter(TotalAbundance > 0)
write.csv(report_neg_marsh, "output/report_neg_marsh.csv")

#greenh
ps_neg_prevalence_greenh <-
  phyloseq::prune_samples("BPCR", ps_greenh) %>%
  plot_abundances()
report_neg_greenh <-
 ps_neg_prevalence_greenh[[1]] %>%
  dplyr::filter(TotalAbundance > 0)
write.csv(report_neg_greenh, "output/report_neg_greenh.csv")
