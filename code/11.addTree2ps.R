###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# August 2021
###.............................................................................
#GOAL: add tree to filtered phyloseq objects
#PROJECT: marsh_metabarcoding
###.............................................................................
library(tidyverse)
library(phyloseq)
library(phangorn)
library(metagMisc)

#read phyloseq objects
ps <-
  readRDS("data/intermediate/ps_filt.rds")

# add trees to phyloseq objects
# read tree
ps_t <-
  phyloseq::read_tree("data/intermediate/raxml/RAxML_bipartitions.marsh.tree") %>%
  phangorn::midpoint() %>%
  phyloseq::merge_phyloseq(ps, .)

ps_t_noRep <-
  phyloseq::subset_samples(ps_t, rep < 2) %>%
  metagMisc::phyloseq_filter_taxa_rel_abund(0) #rm taxa with 0 abundance

# save phyloseq with tree
saveRDS(ps_t, "data/intermediate/ps_t.rds")
saveRDS(ps_t_noRep, "data/intermediate/ps_t_noRep.rds")
