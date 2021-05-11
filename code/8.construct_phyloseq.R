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
taxid_marsh <- readRDS("data/intermediate/taxid_marsh.rds")
taxid_greenh <- readRDS("data/intermediate/taxid_greenh.rds")

# read OTU tables: loads
load("data/intermediate/seqtabNoC.Rdata")

#1.1. read metadata
meta <-
  read.csv("data/raw/metadata.csv", sep = ";")

#assert all samples in the OTU table are present in metadata
assertthat::assert_that(
  all((rownames(seqtabNoC_greenh) %in% meta$sample_name) &
      (rownames(seqtabNoC_marsh) %in% meta$sample_name)))

#2. construct phyloseq object
#phyloseq allows integrating taxonomy table, ASVs table, metadata
#phylogenetic trees into an unique phyloseq object for further manipulation.
# marsh
meta_marsh <-
  dplyr::filter(meta, dataset == "marsh" | sample_name == "BPCR")
rownames(meta_marsh) <- meta_marsh$sample_name

ps_marsh <-
  phyloseq::phyloseq(
    otu_table(seqtabNoC_marsh, taxa_are_rows = FALSE),
    sample_data(meta_marsh),
    tax_table(taxid_marsh)) %>%
    # copy ASV DNA seq to refseq slot and simplify otu_table names
    taxanames2refseq_ps()

# greenhouse
meta_greenh <-
  dplyr::filter(meta, dataset == "greenhouse" | sample_name == "BPCR")
rownames(meta_greenh) <- meta_greenh$sample_name

ps_greenh <-
  phyloseq::phyloseq(
    otu_table(seqtabNoC_greenh, taxa_are_rows = FALSE),
    sample_data(meta_greenh),
    tax_table(taxid_greenh)) %>%
    # copy ASV DNA seq to refseq slot and simplify otu_table names
    taxanames2refseq_ps()

#save phyloseq objects
saveRDS(ps_greenh, "data/intermediate/ps_greenh.rds")
saveRDS(ps_marsh, "data/intermediate/ps_marsh.rds")
