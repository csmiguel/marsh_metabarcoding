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

library(dplyr)
library(DESeq2)
library(phyloseq)

ps <- readRDS("data/intermediate/ps_filt.rds")

#vector with variables to determine differential abundances
variables <- c("species", "season", "rhizosphere")
assertthat::assert_that(all(variables %in% names(sample_data(ps))))

#format taxa
tax <-
  tax_table(ps) %>%
  data.frame() %>%
  tibble::rownames_to_column("seq")

#format OTUs
otut <-
  t(otu_table(ps)) %>%
  data.frame() %>%
  tibble::rownames_to_column("seq")

# geometric means prior to estimate size factors
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# infere differential taxa
deseq <-
  variables %>%
  lapply(function(varx) {
    meta_deseq <-
      phyloseq::phyloseq_to_deseq2(ps,
                  as.formula(paste("~ ", varx, sep = "")))
    geoMeans <- apply(counts(meta_deseq), 1, gm_mean)
    meta_deseq <-
      estimateSizeFactors(meta_deseq, geoMeans = geoMeans) %>%
      DESeq2::DESeq()
    p1 <-
      plotDispEsts(meta_deseq)
  sigtab <-
    as(results(meta_deseq), "data.frame")
  sigtab <- cbind("seq" = rownames(sigtab), sigtab, stringsAsFactors = F)

  #ASV with most contribute to differences among the given condition
  taxa_diff <-
    sigtab[which(sigtab$padj < 0.05),] %>%
    arrange(padj)
  #join taxonomic assignations + frequencies
    merged_results <-
      left_join(taxa_diff, tax, by = "seq") %>%
      left_join(otut, by = "seq")
    return(merged_results)
    }) %>% setNames(variables)

#save differential taxa
saveRDS(deseq, file = "data/intermediate/deseq2.rds")
