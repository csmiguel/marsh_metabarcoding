
plot_deseq2 <- function(ds_result = NULL, out_file = NULL, p_thrs = 0.05,
                        clustr = 1, clustc = 1, use_id = F) {
  # ds_result, is a list of 1 or multiple DESeqResults objects
  # out_file, is the filepath to write the plot in pdf format
  # var2test, is the variable tested in the
  require(plyr); require(DESeq2); require(tidyverse)
  require(pheatmap); require(gtools)
  assertthat::assert_that(
    is.list(ds_result) &
      all(sapply(ds_result, class) == "DESeqResults"),
  msg = "ds_result, is not a list or does not contain DESeqResults objects")

  dsi <-
    plyr::ldply(ds_result, function(x) {
      hdf <-
        as.data.frame(x) %>%
        dplyr::select(log2FoldChange, padj) %>%
        tibble::rownames_to_column("asv") %>%
        dplyr::mutate(contrast = x@elementMetadata@listData$description[4] %>%
                        stringr::str_remove("Wald statistic: [aA-zZ]++ "))
    })
    dsi <-
      if(use_id) mutate(dsi, contrast = paste(.id, contrast, sep = ":")) else dsi
    dsi <-
    dplyr::left_join(dsi, tax, by = c("asv" = "OTU")) %>% # join taxonomy
    {
      log2df <- # data frame with log2 values, wide format
        reshape2::dcast(., phyfam ~ contrast, value.var = "log2FoldChange")
      padjdf <- # data frame with p values, wide format
        reshape2::dcast(., phyfam ~ contrast, value.var = "padj")
      asvs2retain <- # T/F vector with families to retain: with any p value < thr
        select(padjdf, -phyfam) %>%
        apply(1, function(x) !all(x > p_thrs)) %>% {.[is.na(.)] <- F; .}
      log2df <- log2df[asvs2retain, ] # filter wide formated df with log2 abund
      # wide formatted p values df
      padjdf <<- padjdf[asvs2retain, ]  %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("phyfam") %>%
        round(2)# filter
      log2df # ouput log2df to pipe
    } %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("phyfam")
    #replace p values with stars
    padjdf <- gtools::stars.pval(as.matrix(padjdf))

    pheatmap::pheatmap(dsi,
                       fontsize = 5,
                       display_numbers = padjdf,
                       number_format = "%.1f",
                       na_col = "grey",
                       border_color = "grey",
                       cellwidth = 15, cellheight = 6,
                       cutree_rows = clustr,
                       cutree_cols = clustc,
                       labels_col = names(dsi),
                       angle_col = 90,
                       filename = out_file)
    }
