alpha_diversity <- function(ps = NULL, outplot = NULL) {
  require(phyloseq)
  #ps, is a phyloseq object
  #outplot, is the path to a plot for assessing the correlation between
  # no of reads and diversity.
  set.seed(123)
  ps_rarefied <-
    phyloseq::rarefy_even_depth(physeq = ps,
                                rngseed = 1,
                                sample.size = 0.9 * min(sample_sums(ps)))

  #estimate alpha diversity
  alpha_div <-
    phyloseq::estimate_richness(ps_rarefied) %>%
    dplyr::mutate(Pielou = Shannon/log(Observed)) %>% #compute Pielou diversity
    `rownames<-`(sample_names(ps_rarefied)) # fix rownames
  #check outliers
  pdf(outplot)
  plot(phyloseq::sample_sums(ps), alpha_div$Shannon,
    xlab = "Number of reads",
    ylab = "Shannon H'",
    type = "n")
  text(phyloseq::sample_sums(ps), alpha_div$Shannon,
       labels = phyloseq::sample_names(ps_rarefied),
       cex = 0.5, col = "blue")
  dev.off()
  return(alpha_div)
}
