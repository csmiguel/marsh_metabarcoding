###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: create dna alignments for tree reconstruction
#PROJECT: marsh_metabarcoding
###.............................................................................
library(DECIPHER)
library(phyloseq)
library(dplyr)
library(Biostrings)
library(phangorn)

#phyloseq
ps <- readRDS("data/intermediate/ps_filt.rds")

# create dir raxml
dirraxml <- "data/intermediate/raxml"
if(!dir.exists(dirraxml)) dir.create(dirraxml)

# align sequences, and write to phylip
phyloseq::refseq(ps) %>%
  DECIPHER::AlignSeqs(anchor = NA) %>%
  Biostrings::DNAMultipleAlignment() %>%
  phangorn::as.phyDat() %>%
  phangorn::write.phyDat(file = file.path(dirraxml, "marsh.phy"),
                         format = "phylip")
