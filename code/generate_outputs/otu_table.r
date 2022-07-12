###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: write otu table to csv and phyloseq Rdata
#PROJECT: marsh_metabarcoding
###.............................................................................
#
library(phyloseq)
library(tidyverse)
library(tibble)

# read filtered phyloseq
ps <-
  readRDS("data/intermediate/ps_t_noRep.rds")

# write phyloseq to otu table
# create data frame for export
otu <-
  as.data.frame(t(otu_table(ps))) %>%
  tibble::rownames_to_column("asv")

tax <-
  as.data.frame(tax_table(ps)) %>%
  tibble::rownames_to_column("asv")

refsq <-
  refseq(ps) %>%
  as.data.frame() %>%
  rename(dna_sequence = x) %>%
  tibble::rownames_to_column("asv")

tax_otu_refseq <-
  dplyr::left_join(tax, refsq) %>%
  dplyr::left_join(otu)

# save data frame as csv
write.csv(tax_otu_refseq,
          file = "output/otu_table.csv",
          quote = T)

# save phyloseq to output for data availability
save(ps, file = "output/phyloseq.Rdata")
