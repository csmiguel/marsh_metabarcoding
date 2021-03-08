###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# March 2021
###.............................................................................
#GOAL: Determination of ASVs
#We followed the strategy described in [Callahan et al. 2016](http://dx.doi.org/10.12688/f1000research.8986.2)
# in its latest [update](https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html)
#PROJECT: spartina-metarizo
###.............................................................................
library(dplyr)
library(dada2)
#Set path to path directory with trimmed sequences from cutadapt
path <- "data/intermediate"

#select samples from the field: they start with "S".
fnFs <-
  sort(list.files(path,
    pattern = "S.*_R1.cutadapt.fastq.gz", full.names = TRUE))
fnRs <-
  sort(list.files(path,
    pattern = "S.*_R2.cutadapt.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#create path for sequences that will be filtered downstream
filt_path <- file.path(path, "filtered")
if (!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#save objects
# path to trimmed reads
saveRDS(fnFs, "data/intermediate/fnFs.rds")
saveRDS(fnRs, "data/intermediate/fnRs.rds")
