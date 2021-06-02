###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: create dna alignments for tree reconstruction
#PROJECT: spartina-metarizo
###.............................................................................
library(DECIPHER)
library(phyloseq)
library(dplyr)
library(Biostrings)
library(phangorn)

#marsh
ps_marsh <- readRDS("data/intermediate/ps_marsh_filt.rds")
# create dir raxml
dir_marsh <- "data/intermediate/raxml_marsh"
if(!dir.exists(dir_marsh)) dir.create(dir_marsh)
# align sequences, and write to phylip
phyloseq::refseq(ps_marsh) %>%
  DECIPHER::AlignSeqs(anchor = NA) %>%
  Biostrings::DNAMultipleAlignment() %>%
  phangorn::as.phyDat() %>%
  phangorn::write.phyDat(file = file.path(dir_marsh, "marsh.phy"),
                         format = "phylip")
#greenh
ps_greenh <- readRDS("data/intermediate/ps_greenh_filt.rds")
# create dir raxml
dir_greenh <- "data/intermediate/raxml_greenh"
if(!dir.exists(dir_greenh)) dir.create(dir_greenh)
# align sequences, and write to phylip
phyloseq::refseq(ps_greenh) %>%
  DECIPHER::AlignSeqs(anchor = NA) %>%
  Biostrings::DNAMultipleAlignment() %>%
  phangorn::as.phyDat() %>%
  phangorn::write.phyDat(file = file.path(dir_greenh, "greenh.phy"),
                         format = "phylip")
