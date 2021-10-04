deseqp <- function(res, pval) {
  # res is an object from DESeq2::results
  # pval is the threshold of pvalues for filtering asvs
  require(dplyr)
    h <- as(res, "data.frame")
  #
    sigtab <-
      cbind("seq" = rownames(h),
            h, stringsAsFactors = F)

    #ASV with most contribute to differences among the given condition
    sigtab[which(sigtab$padj < pval), ] %>%
      arrange(padj)
      }
