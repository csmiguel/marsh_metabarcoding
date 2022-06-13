###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# Feb 2022
###.............................................................................
#GOAL: assign taxonomy to ASVs using Silva v138 database
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(ShortRead)
library(dplyr)
library(tibble)

# read phyloseq
ps <- readRDS("data/intermediate/ps_t_noRep.rds")

# extract sequences
seqs <- phyloseq::refseq(ps)

# 1. write sequences to fasta file
ShortRead::writeFasta(seqs,"genbank-submission/asvs.fa")
#>ASV3234
#TGGGGAATATTGCACAATGGAGGAAACTCTGATGCAGCAACGCCGCGTGGAGGATGACACTTTTCGGAGCGTAAACTCCT
#...

# 2. create BioSample mapping table required in the submission wizard.

#   GenBank requires the following mapping file. As indicated in an email from GenBank admin I must create a file
#   with 2 columns, and only 1 BioSample per ASV. However, each ASV can be present in multiple BioSamples. The
#   solution that GenBank gave me is to submit the mapping file through the wizard with only one BioSample per
#   ASV and then another file by email to gb-admin@ncbi.nlm.nih.gov with all the BioSample separated by commas.
# File 1: 
# Sequence_ID	biosample_accession
# ASV3234 BioSampleID1
# ASV6728 BioSampleID2
#
# File 2:
# Seq1	BioSampleID1,BioSampleID2
# Seq2	BioSampleID1,BioSampleID3
#
#   2.1 gather biosample info:
meta <-
  read.csv2("genbank-submission/metadata-9191966-processed-ok.tsv",
            sep = "\t") %>%
  dplyr::select(sample_name, biosample_accession)

# data frame otu table
otut <-
  t(as.data.frame(otu_table(ps)@.Data))
# replace sample names by BioSample accession
colnames(otut) <- meta$biosample_accession[match(colnames(otut), meta$sample_name)]

s_names <- colnames(otut) # sample names

# File 1: all ASVs and the sample name of their first occurrence.
apply(otut, 1, function(x) {
  # return first sample in which the ASV is present
  s_names[match(T, x > 0)]
  }
) %>%
  as.data.frame %>%
  tibble::rownames_to_column("Sequence_ID") %>%
  dplyr::rename("biosample_accession" = 2) %>%
  write.table("genbank-submission/biosamples_wizard.txt",
              sep = "\t",
              quote = F,
              row.names = F)

# File 2: all ASVs and all BioSamples.
apply(otut, 1, function(x) {
  # return all sample in which the ASV is present
  s_names[x > 0] %>% paste(collapse = ",")
}
) %>%
  as.data.frame %>%
  tibble::rownames_to_column("Sequence_ID") %>%
  dplyr::rename("biosample_accession" = 2) %>%
  write.table("genbank-submission/biosamples_email.txt",
              sep = "\t",
              quote = F,
              row.names = F)
