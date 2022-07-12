###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: truncate reads
#PROJECT: marsh_metabarcoding
###.............................................................................
library(dplyr)
library(dada2)
#load filtering parameters
source("code/parameters/trim.r")

#path to F and R trimmed reads
fnFs <- readRDS("data/intermediate/fnFs.rds")
fnRs <- readRDS("data/intermediate/fnRs.rds")

#vector with sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Set path to path directory with trimmed sequences from cutadapt
path <- "data/intermediate"

#create path for sequences that will be filtered downstream
filt_path <- file.path(path, "filtered")
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Based on the quality plots we decided to truncate R1 reads to 270 nt and R2
#reads to 185 nt. Considering the amplicon length is ~425 nt (without primers),
#we would still have have ~30 nt of overlap for merging R1 and R2, downstream.

assertthat::assert_that(all(
  ((trunc_f + trunc_r) - amplicon_length) / 2 > 12, # at least this overlap
  trunc_f < read_length, # trunc length is below read length
  trunc_r < read_length),
  msg = "check read length and truncation paramters")

truncated <-
  dada2::filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    maxN = 0,
    maxEE = expected_errors,
    truncQ = 2,
    compress = TRUE,
    multithread = TRUE,
    truncLen = c(trunc_f, trunc_r))

#save objects
write(sample.names, "data/intermediate/sample-names")
saveRDS(truncated, "data/intermediate/truncated_reads.rds")
saveRDS(sample.names, "data/intermediate/sample_names.rds")
save(filt_path, filtFs, filtRs, file = "data/intermediate/filt.Rdata")
