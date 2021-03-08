###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# March 2021
###.............................................................................
#GOAL: generate quality reports from trimmed reads.
#PROJECT: spartina-metarizo
###.............................................................................

library(fastqcr)

fnFs <- readRDS("data/intermediate/fnFs.rds")
fnRs <- readRDS("data/intermediate/fnRs.rds")

#create fastqc resports
fastqcr::fastqc(fq.dir = "data/intermediate/",
  qc.dir = "data/intermediate/fastqc",
  threads = 2,
  fastqc.path = "/usr/local/bin/fastqc")

#create quality plots
qualplt <-
  dada2::plotQualityProfile(c(fnFs, fnRs),
  n = 1000, aggregate = F)

sizef <- ceiling(sqrt(length(c(fnFs, fnRs))))
ggsave("output/quality-plots.pdf", plot = qualplt, device = "pdf",
  width = 7 * sizef,
  height = 5 * sizef,
  units = "cm")
