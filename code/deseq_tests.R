###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: infere taxa with differential abundances by selected variables
#PROJECT: spartina-metarizo
###.............................................................................
library(tidyverse)
library(DESeq2)
library(phyloseq)

# read phyloseq object and remove replicas
ps <-
  readRDS("data/intermediate/ps_t.rds") %>%
  phyloseq::subset_samples(rep < 2)
# load funcions
source("code/functions/gm_means.r")

#####TESTS#########
# I will do multiple tests:
# 1. effect of plant species on the microbial community: community ~ species
# 2. per species and across stations, are there differences between rhizosphere
# and naked soil: community ~ species:rhizosphere
# 3.effect of season on species: community ~ species:season
##################

# first step is to transform the phyloseq object into a deseq2 input
# convert phyloseq 2 deseq with formula. All variables in the formula are
# converted to factors. Then, it is possible to change the formula with "design"
temp1 <-
      phyloseq::phyloseq_to_deseq2(ps,
        as.formula("~ species + season + rhizosphere"))
# compute geometric means for the counts
geoMeans <- apply(counts(temp1), 1, gm_mean)

# deseq2 analysis. By estimating size factors each ASV and sample abundances are
# weighted for computing the statistics
ds <- estimateSizeFactors(temp1, geoMeans = geoMeans)

# TEST 1
temp2 <- ds
design(temp2) <- ~species
ds_test1 <- DESeq2::DESeq(temp2)



#Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level. The comparison is printed at the top of the output: dex trt vs untrtå
sigtab <-
  DESeq2::results(meta_deseq) %>%
      as("data.frame")

    sigtab <- cbind("seq" = rownames(sigtab), sigtab, stringsAsFactors = F)

    #ASV with most contribute to differences among the given condition
    taxa_diff <-
      sigtab[which(sigtab$padj < 0.05),] %>%
      arrange(padj)
