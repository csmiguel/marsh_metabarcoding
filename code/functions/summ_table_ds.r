# function to create a summary table from deseq results
summ_table_ds <- function(ds = ds, pthres = 1) {
  # ds are tidy deseq restuls
  # pthres is the threshold on the adjusted pvalues to include in the table
  ds %>%
    dplyr::filter(padj < pthres) %>%
    group_by(test, contrast) %>%
    summarise(number = n(),
              pos = sum(log2FoldChange > 0),
              neg = sum(log2FoldChange < 0)) %>%
    dplyr::mutate(factor = plyr::mapvalues(test,
                                           c("test1", "test2", "test3", "test4"),
                                           c("species_plot", "rhizosphere", "rhizosphere", "season")),
                  contrast = stringr::str_remove(contrast, "contrast [1-9][0]? ")) %>%
    tidyr::separate(contrast, c("level_1", "level_2"), sep = " vs ") %>%
    tidyr::separate(level_1, c("strata", "level_1"), " ", fill = "left") %>%
    as.data.frame %>%
    dplyr::select(factor, strata, level_1, level_2, number, pos, neg) %>%
    dplyr::rename(total = number)
}