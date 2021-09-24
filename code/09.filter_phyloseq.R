###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: construct phyloseq objects
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(dplyr)

#read phyloseq objects
ps <- readRDS("data/intermediate/ps.rds")

#load functions for filtering phyloseq objects
dir("code/functions", "ps_filter", full.names = T) %>% sapply(source)

#open connection
sink(file = "output/filter_phyloseq.txt")

# filter ASVs from phyloseq
ps_filt <-
  ps_filter_phylum_is_NA(ps) %>%
    ps_filter_contaminants("BPCR") %>%
    ps_filter_organelles() %>%
    ps_filter_prevalence(
      mult_threshold = log(nsamples(ps)) * 50, #I use log to account for
      # the fact that with larger sample size the chances of having shared false
      # positivies across samples could increase in a logarithmic manner.
      prevalence_threshold = 2
                        ) %>%
    ps_filter_relative_abundance(mean_prop = 5e-5)  #I use a more relaxed
    #threshold for marsh plants because since they come
    # from different species I expect them to have lower mean_prop in global
sink()

#save filtered phyloseq
saveRDS(ps_filt, "data/intermediate/ps_filt.rds")
