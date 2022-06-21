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
