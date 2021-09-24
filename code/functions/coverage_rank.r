# function to compute the coverage per taxnomic rank. Coverage is defined here as the number of reads assigned to a given non-NA taxonomic rank divided by the total number of reads.
rank_cov <- function(phyloseq_object) {
  require(rlang)
  require(phyloseq)
  require(dplyr)
  
  #melt phyloseq object
  psmelted <- phyloseq::psmelt(phyloseq_object)
  
  #total no of reads
  total_reads <- sum(psmelted$Abundance)
  
  # rank names except species
  rn <- rank_names(ps)[!grepl("species", rank_names(ps))]
  
  # compute proportion of non_NA tax ranks
  rn %>%
    sapply(function(x) {
      psmelted %>%
        dplyr::filter(!is.na(!!!rlang::parse_exprs(x))) %>%
        pull(Abundance) %>%
        sum
    }) %>%
    {. / total_reads} %>% round(2)
  }