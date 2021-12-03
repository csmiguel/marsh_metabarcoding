# generate outputs:
# Venn diagram for comparisson between rhizosphere vs  non-rhizospheric soil.
# table with number of differential families associated to each condition
# table with pairwise comparisson of fammilies between species plot.
library(tidyverse)
library(ggVennDiagram)
library(gridExtra)
library(flextable)

# load function for Venn diagrams
source("code/functions/venn_ds.r")
source("code/functions/summ_table_ds.r")

# load deseq tidy data
ds <-
  readRDS("data/intermediate/deseq_results_merged.rds") %>%
  dplyr::mutate(contrast = gsub(pattern = "Salicornia ramosissima", replacement = "Salicornia_ramosissima", contrast) %>%
                  gsub(pattern = "rocaulon macr", replacement = "rocaulon_macr") %>%
                  gsub(pattern = "Spartina ", replacement = "Spartina_")) # homogeneize names of plants

#1. summary table
# summary table for significant families
sum_05 <- summ_table_ds(ds, pthres = 0.05)

# summary table for all data (this is done so as to show as well those levels which do not have sigficant families)
sum_all <- summ_table_ds(ds, pthres = 1)

# number of differential groups per test
diff_test <-
  ds %>%
  dplyr::filter(padj < .05) %>%
  dplyr::select(test, family) %>%
  dplyr::mutate(factor = gsub("test2", "test3", test)) %>%
  distinct() %>%
  group_by(factor) %>%
  summarise(total_per_factor = n()) %>%
  dplyr::mutate(factor = plyr::mapvalues(factor,
               c("test1", "test2", "test3", "test4"),
               c("species_plot", "rhizosphere", "rhizosphere", "season")))

# merge into one table
summ_table <-
  dplyr::left_join(dplyr::select(sum_all, -total, -pos, -neg),
                   sum_05) %>%
  dplyr::mutate(total = tidyr::replace_na(total, 0),
                pos = tidyr::replace_na(pos, 0),
                neg = tidyr::replace_na(neg, 0)) %>%
  dplyr::left_join(diff_test)

# format flextable
flextable(summ_table) %>%
  flextable::save_as_docx(path = "output/summ_table_ds.docx")

# 2. Venn plot for rhizospheric soil vs bulk soil. The goal is to visualize if there are specific families associated to rhizoshpere.
# plot negative values
p_test3_neg <- ds_ven(ds = ds,
                      testi = "test3",
                      logchange_pos = F,
                      col_high = "grey")
# plot positive values
p_test3_pos <- ds_ven(ds = ds,
                      testi = "test3",
                      logchange_pos = T,
                      col_high = "grey")

# save plot
ggsave(grid.arrange(p_test3_pos, p_test3_neg),
       width = 5,
       height = 7,
       file = "output/venn_rhizosphere.pdf")

# 3. diff table species plot (not included in manuscript)
ds_sp_pos <-
  ds %>%
    dplyr::filter(test == "test1" & padj < 0.05 & log2FoldChange > 0) %>%
    tidyr::separate(col = contrast, into = letters[1:5], sep = " ") %>%
    dplyr::rename(sp1 = c, sp2 = e) %>%
    dplyr::select(sp1, sp2, family) %>%
    reshape2::dcast(formula = sp1 ~ sp2)

ds_sp_neg <-
  ds %>%
  dplyr::filter(test == "test1" & padj < 0.05 & log2FoldChange < 0) %>%
  tidyr::separate(col = contrast, into = letters[1:5], sep = " ") %>%
  dplyr::rename(sp1 = c, sp2 = e) %>%
  dplyr::select(sp1, sp2, family) %>%
  reshape2::dcast(formula = sp1 ~ sp2)
    
xx <-   ds %>%
  dplyr::filter(test == "test1" & padj < 0.05 & log2FoldChange < 0) %>%
  tidyr::separate(col = contrast, into = letters[1:5], sep = " ") %>%
  dplyr::rename(sp1 = c, sp2 = e) %>%
  dplyr::select(sp1, sp2, family) %>%
  group_by(sp1, sp2) %>% summarise(total = n())
