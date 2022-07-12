###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# August 2021
###.............................................................................
#GOAL: FAPROTAX input creation
#PROJECT: marsh_metabarcoding
###...........................................................................
library(phyloseq)
library(dplyr)
library(reshape2)

#read phyloseq
ps <- readRDS("data/intermediate/ps_t_noRep.rds")
h <- psmelt(ps)

hh <- dcast(h, OTU + phylum + class + order + family + genus ~ sample_name,
      value.var = "Abundance")
df <-
  hh %>%
  mutate(taxonomy = paste(phylum, class, order, family, genus, sep = "; ") %>%
           gsub(pattern = "; NA*", replacement = "")) %>%
  select(-phylum, -order, -class, -family, -genus, -OTU)
write.table(df, file = "data/intermediate/faprotax.tsv", quote = F, sep = "\t", row.names = F)
