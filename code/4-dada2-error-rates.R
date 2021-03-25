###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# March 2021
###.............................................................................
#GOAL: estimate error rates
# estimate error rates of F and R from subsample sequences
#PROJECT: spartina-metarizo
###.............................................................................
# the computing power needed growths quadratically with # of seqs
# [big data](https://benjjneb.github.io/dada2/bigdata.html).
# For that reason, we calculated error matrices from random subsamples.
library(dplyr)
library(ShortRead)
library(dada2)

#load data
sample.names <- readRDS("data/intermediate/sample_names.rds")
load("data/intermediate/filt.Rdata") #load filtFs, filtRs and filt_path
n <- 10000 # number of reads to subsample to estimate the error rates.

#sampler
seq_along(sample.names) %>%
  lapply(function(x) {
    sample_name <- sample.names[x] #sample name
    #path to files
    fqf <- grep(pattern = paste0(sample.names[x], "_"), x = filtFs, value = T)
    fqr <- grep(pattern = paste0(sample.names[x], "_"), x = filtRs, value = T)
    assertthat::assert_that(length(fqf) == 1 & length(fqr) == 1)
    #sample
    #forward
    subsamplef <- #subsample reads
      ShortRead::FastqSampler(fqf, n) %>%
      ShortRead::yield(sampler)

    pathf <- file.path(filt_path,
                       paste0(sample_name,
                       "_F_filt_10000_fastq.gz"))
    ShortRead::writeFastq(subsamplef, pathf)
    #reverse
    subsampler <- #subsample reads
      ShortRead::FastqSampler(fqr, n) %>%
      ShortRead::yield(sampler)

    pathr <- file.path(filt_path,
                       paste0(sample_name,
                       "_R_filt_10000_fastq.gz"))
    ShortRead::writeFastq(subsampler, pathr)
    })

#lear error rates
errF <-
  list.files(filt_path, pattern = "F_filt_10000_fastq.gz", full.names = T) %>%
  dada2::learnErrors(multithread = TRUE)
errR <-
  list.files(filt_path, pattern = "R_filt_10000_fastq.gz", full.names = T) %>%
  dada2::learnErrors(multithread = TRUE)

#plot errors
plot_errF <- dada2::plotErrors(errF, nominalQ = TRUE)
plot_errR <- dada2::plotErrors(errF, nominalQ = TRUE)

ggplot2::ggsave("output/plot-errorsF.pdf", plot_errF)
ggplot2::ggsave("output/plot-errorsR.pdf", plot_errR)

#save errors
save(errF, errR, file = "data/intermediate/errors.Rdata")
