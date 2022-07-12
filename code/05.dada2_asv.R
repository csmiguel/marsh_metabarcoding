###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: dada ASV determination
#PROJECT: marsh_metabarcoding
###.............................................................................
library(dada2)
library(dplyr)
load("data/intermediate/errors.Rdata")
load("data/intermediate/filt.Rdata")

#de-replication
# removes identical sequences, lowering the computing burden downstream
derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)

#rename dereps
names(derepFs) <-
  names(derepFs) %>%
    stringr::str_remove("_F_filt.fastq.gz")

names(derepRs) <-
  names(derepRs) %>%
    stringr::str_remove("_R_filt.fastq.gz")

#estimate ASVs
source("code/functions/dada-custom-functions.r")

dadaFs <- dada_1(derepFs, "S-[0-9]|BPCR", errF)
dadaRs <- dada_1(derepRs, "S-[0-9]|BPCR", errR)

#save objects
save(dadaFs, dadaRs,
  file = "data/intermediate/dada.Rdata")
save(derepFs, derepRs,
  file = "data/intermediate/dereps.Rdata")
