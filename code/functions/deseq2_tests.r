ds_test <- function(ps = NULL, tests = NULL) {
 # ps, phyloseq object
 # tests, a list of contrasts. Each contrast is a vector where [1] is the
 #       formula, [2] is the numerator for deseq2 results and [3] is the denominator.

  lapply(tests,
    function(x) {
        ds <-
              phyloseq::phyloseq_to_deseq2(ps,
                as.formula(paste0("~", x[1])))
        # compute geometric means for the counts
        geoMeans <- apply(counts(ds), 1, gm_mean)

        ds <-
          estimateSizeFactors(ds, geoMeans = geoMeans) %>%
          DESeq2::DESeq()

        results(ds, contrast = x)
      })
    }
