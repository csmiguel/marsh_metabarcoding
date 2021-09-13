###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: process dada asvs, invoving merge pairs, makeSequenceTable
#   and remove chimeras.
#PROJECT: spartina-metarizo
###.............................................................................
# load data from dada2::dada
load("data/intermediate/dada.Rdata") #warning: high RAM consumption
load("data/intermediate/dereps.Rdata")
# load customized dada2 functions
source("code/functions/dada-custom-functions.r")
source("code/parameters/filtering_ASVs.r")

library(dada2)
library(dplyr)
#merge paired reads. Spurious sequence variants are further
# reduced by merging overlapping reads:
# use custom functions to process both datasets separately
mergers <-
    mergerPairs_1(dada_f = dadaFs,
                  dada_r = dadaRs,
                  derep_f = derepFs,
                  derep_r = derepRs,
                  pattern_sampleset = "S-[0-9]|BPCR")

#construct sequence table (analogous to an OTU table)
seqtab <- dada2::makeSequenceTable(mergers)

#number of samples and variants
dim(seqtab)

#distribution of sequence lengths
cat("Sequence length distribution of ASVs before removing chimeras")
cat("\nFor marsh plants:")
table(nchar(getSequences(seqtab)))
cat("\nWe retained reads in the range", range(filtering_seqtab))

#requences that are much longer or shorter than expected may be the result of
# non-specific priming, and may be worth removing
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% filtering_seqtab]

#remove bimeras
seqtabNoC <- dada2::removeBimeraDenovo(seqtab2)

#save objects
#merged seqs
saveRDS(mergers,
  file = "data/intermediate/mergers.rds")
#merged seqs size filtered
saveRDS(seqtab2,
    file = "data/intermediate/seqtab2.rds")
#OTU tables no-chimeras
saveRDS(seqtabNoC,
     file = "data/intermediate/seqtabNoC.rds")
