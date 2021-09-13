alpha_diversity_bt <- function(ps = NULL, bt = NULL) {
  #ps is a phyloseq object
  #bt is the number of bootstrap replicates
  seq_len(bt) %>%
    sapply(function(x) {
      phyloseq::rarefy_even_depth(physeq = ps,
                                  sample.size = 0.9 * min(sample_sums(ps))) %>%
      phyloseq::estimate_richness(measures = "Shannon")}) %>%
      do.call(what = cbind) %>%
      as.data.frame() %>%
      setNames(paste0("bt", seq_len(bt))) %>%
      {cbind(as(sample_data(ps), "data.frame"), .)}
      }
