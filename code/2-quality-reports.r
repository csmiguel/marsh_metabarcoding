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
library(dada2)
library(dplyr)

path <- "data/intermediate"
#select samples from the field: they start with "S".
fnFs <-
  sort(list.files(path,
    pattern = "R1.cutadapt.fastq.gz", full.names = TRUE))
fnRs <-
  sort(list.files(path,
    pattern = "R2.cutadapt.fastq.gz", full.names = TRUE))

#create fastqc resports on trimmed fastq
# remove potential files stored in tempdir
list.files(tempdir(), pattern = "fastq", full.names = T) %>% file.remove()
# softlink of fastq files to a temporary directory for fastqc
file.symlink(from = file.path(getwd(),
                              list.files(path,
                    pattern = "R[1-2].cutadapt.fastq.gz", full.names = T)),
             to =  tempdir())
# execute fastq
fastqcr::fastqc(fq.dir = tempdir(),
  qc.dir = "output/fastqc",
  threads = 2,
  fastqc.path = "/usr/local/bin/fastqc")

#create quality plots
temp1 <- list(fnFs = fnFs, fnRs = fnRs)
seq_along(temp1) %>%
  lapply(function(x) {
    qualplt <-
      dada2::plotQualityProfile(temp1[[x]],
      n = 1000, aggregate = F)
    # save quality plots
    sizef <- ceiling(sqrt(length(temp1[[x]])))
    ggsave(paste0("output/quality-plots-", names(temp1)[x], ".pdf"),
      plot = qualplt, device = "pdf",
      width = 7 * sizef,
      height = 5 * sizef,
      units = "cm")
    })

# path to trimmed reads
saveRDS(fnFs, "data/intermediate/fnFs.rds")
saveRDS(fnRs, "data/intermediate/fnRs.rds")
