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
ps_marsh <- readRDS("data/intermediate/ps_marsh.rds")
ps_greenh <- readRDS("data/intermediate/ps_greenh.rds")

#load functions for filtering phyloseq objects
dir("code/functions", "ps_filter", full.names = T) %>% sapply(source)
#open connection
sink(file = "output/filter_phyloseq.txt")
# filter ASVs from phyloseq
#marsh
ps_marsh_filt <-
  ps_filter_phylum_is_NA(ps_marsh) %>%
    ps_filter_contaminants("BPCR") %>%
    ps_filter_organelles() %>%
    ps_filter_prevalence(
      mult_threshold = log(nsamples(ps_marsh)) * 50, #I use log to account for
      # the fact that with larger sample size the chances of having shared false
      # positivies across samples could increase in a logarithmic manner.
      prevalence_threshold = 2
                        ) %>%
    ps_filter_relative_abundance(mean_prop = 5e-5)  #I use a more relaxed
    #threshold for marsh plants because since they come
    # from different species I expect them to have lower mean_prop in global

#greenhouse
ps_greenh_filt <-
  ps_filter_phylum_is_NA(ps_greenh) %>%
      ps_filter_contaminants("BPCR") %>%
      ps_filter_organelles() %>%
      ps_filter_prevalence(
          mult_threshold = log(nsamples(ps_greenh)) * 50,
          prevalence_threshold = 2
      ) %>%
      ps_filter_relative_abundance(mean_prop = 1e-4)
#write to file and close connection
sink()
#save filtered phyloseq
saveRDS(ps_marsh_filt, "data/intermediate/ps_marsh_filt.rds")
saveRDS(ps_greenh_filt, "data/intermediate/ps_greenh_filt.rds")
