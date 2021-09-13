###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# May 2021
###.............................................................................
#GOAL: construct phyloseq objects
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(dplyr)

source("code/functions/taxanames2refseq_ps.r")

#1. read data
# read seqs with taxonomy
taxid <- readRDS("data/intermediate/taxid.rds")

# read OTU tables: loads
seqtabNoC <- readRDS("data/intermediate/seqtabNoC.rds")

#1.1. read metadata
meta <-
  read.csv("data/raw/metadata.csv", sep = ";")

#assert all samples in the OTU table are present in metadata
assertthat::assert_that(
  all((rownames(seqtabNoC) %in% meta$sample_name)))

#2. construct phyloseq object
#phyloseq allows integrating taxonomy table, ASVs table, metadata
#phylogenetic trees into an unique phyloseq object for further manipulation.
rownames(meta) <- meta$sample_name

ps <-
  phyloseq::phyloseq(
    otu_table(seqtabNoC, taxa_are_rows = FALSE),
    sample_data(meta),
    tax_table(taxid)) %>%
    # copy ASV DNA seq to refseq slot and simplify otu_table names
    taxanames2refseq_ps()

#save phyloseq objects
saveRDS(ps, "data/intermediate/ps.rds")
