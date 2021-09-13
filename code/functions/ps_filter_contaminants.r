ps_filter_contaminants <- function(ps = NULL, blanks = NULL) {
  #ps is a phyloseq object
  #blanks is a character vector with sample name/s of blanks
  require(decontam)
  require(phyloseq)
  #assertions
  assertthat::assert_that(all(blanks %in% sample_names(ps)))
  #blanks
  blank_postions <- rownames(otu_table(ps)) %in% blanks
  #contanimanion analysis
  contam_df <- decontam::isContaminant(ps, neg = blank_postions)
  # contaminant ASVs
  contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
  # non-contaminant ASVs
  non_contam_asvs <- row.names(contam_df[contam_df$contaminant == FALSE, ])
  h <-
    prune_taxa(non_contam_asvs, ps) %>%
    {prune_samples(! sample_names(.) %in% blanks, .)}
  if (identical(contam_asvs, character(0))) {
    cat("\nNo contaminant ASVs were identified.\n")
  } else if (length(contam_asvs) > 0) {
    cat("\n", contam_asvs,
      "was/were identified as contaminants and removed from the phyloseq object.\n")
    }
    cat("\nSample/s",
      sample_names(h), "were kept.\nSample/s", blanks, "was/were removed.\n")
  #report results from contamint filtering
  no_reads <-
    mean(rowSums(otu_table(h)) /
    rowSums(otu_table(ps)[!sample_names(ps) %in% blanks]))
  no_taxa <- ntaxa(h) / ntaxa(ps)
  #report some results from the filtering
  cat("\nAfter removing contaminants and blanks a",
  no_reads, "of the reads and a",
  no_taxa, "of the ASVs were retained.\n")
    return(h)
  }
