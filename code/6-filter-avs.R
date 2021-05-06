###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2021
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
mergers_marsh <-
    mergerPairs_1(dada_f = marsh_dadaFs,
                  dada_r = marsh_dadaRs,
                  derep_f = derepFs,
                  derep_r = derepRs,
                  pattern_sampleset = "S-[0-9]|BPCR")

mergers_greenh <-
    mergerPairs_1(dada_f = greenh_dadaFs,
                  dada_r = greenh_dadaRs,
                  derep_f = derepFs,
                  derep_r = derepRs,
                  pattern_sampleset = "MA-[0-9]|BPCR")

#construct sequence table (analogous to an OTU table)
seqtab_marsh <- dada2::makeSequenceTable(mergers_marsh)
seqtab_greenh <- dada2::makeSequenceTable(mergers_greenh)

#number of samples and variants
dim(seqtab_marsh)
dim(seqtab_greenh)

#distribution of sequence lengths
cat("Sequence length distribution of ASVs before removing chimeras")
cat("\nFor marsh plants:")
table(nchar(getSequences(seqtab_marsh)))
cat("\nFor Spartina maritima on pots:")
table(nchar(getSequences(seqtab_greenh)))
cat("\nWe retained reads in the range", range(filtering_seqtab))

#requences that are much longer or shorter than expected may be the result of
# non-specific priming, and may be worth removing
seqtab_marsh_2 <- seqtab_marsh[, nchar(colnames(seqtab_marsh)) %in% filtering_seqtab]
seqtab_greenh_2 <- seqtab_greenh[, nchar(colnames(seqtab_greenh)) %in% filtering_seqtab]

#remove bimeras
seqtabNoC_marsh <- dada2::removeBimeraDenovo(seqtab_marsh_2)
seqtabNoC_greenh <- dada2::removeBimeraDenovo(seqtab_greenh_2)

#save objects
#merged seqs
save(mergers_marsh, mergers_greenh,
    file = "data/intermediate/mergers.Rdata")
#merged seqs size filtered
save(seqtab_marsh_2, seqtab_greenh_2,
    file = "data/intermediate/seqtab2.Rdata")
#OTU tables no-chimeras
save(seqtabNoC_marsh, seqtabNoC_greenh,
     file = "data/intermediate/seqtabNoC.Rdata")
