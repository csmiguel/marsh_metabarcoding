ds_ven <- function(ds = ds,
                   testi = NULL,
                   padj_filt = .05,
                   logchange_pos = T,
                   col_low = "white",
                   col_high = "red",
                   regex_contrast = NULL) {
  input1 <-
    ds %>%
    dplyr::filter(test == testi & padj < padj_filt &
                    if(logchange_pos) log2FoldChange > 0 else log2FoldChange < 0) %>%
    dplyr::select(contrast, family) %>%
    plyr::dlply(~contrast, function(x) x$family)
  #filter contrasts
  if(!is.null(regex_contrast)) input1 <- input1[grepl(regex_contrast, names(input1))]
  # create Venn data
  venn <- Venn(input1)
  data <- process_data(venn)
  ggVennDiagram(input1,
                label_alpha = 0,
                set_size = 3,
                label_color = "black", edge_color = "black") +
    geom_sf(size = 1, lty = 1, color = 1, data = venn_setedge(data), show.legend = F) +
    scale_fill_gradient(low = col_low, high = col_high)
}
