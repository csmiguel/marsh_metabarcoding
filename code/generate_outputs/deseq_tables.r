# produce deseq tables
library(DESeq2)
library(tidyverse)

# deseq results
ds <- readRDS("data/intermediate/deseq2_results.rds")

# tax table by family
tax <- readRDS("data/intermediate/tax_df.rds")

# merge all results into a unique table
xy <- 
  lapply(ds, function(y) {
    # add names to list in case they are null
    yy <- if(is.null(names(y))) paste("contrast", seq_along(y)) else names(y)
    names(y) <- yy
    plyr::ldply(y, function(x) {
        as.data.frame(x) %>%
        tibble::rownames_to_column("asv") %>%
        dplyr::mutate(contrast = x@elementMetadata@listData$description[4] %>%
                        stringr::str_remove("Wald statistic: [aA-zZ]++ "))
      }) %>%
      dplyr::rename(contrast_n = .id)
    }) %>%
  dplyr::bind_rows(.id = "test") %>%
  dplyr::mutate(contrast = paste(contrast_n, contrast)) %>%
  dplyr::left_join(tax, by = c("asv" = "OTU")) %>%
  dplyr::select(-contrast_n, -asv) %>%
  as_tibble()

# write table to file
write.csv(xy, file = "output/deseq_results.csv", quote = T)

# write table to rds
saveRDS(xy, "data/intermediate/deseq_results_merged.rds")
