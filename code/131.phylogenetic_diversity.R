library(phyloseq)
library(metagMisc)
library(tidyverse)
library(ggplot2)

# read data
ps <- readRDS("data/intermediate/ps_t_noRep.rds")

set.seed(123)
ps_rarefied <-
  phyloseq::rarefy_even_depth(physeq = ps,
                              rngseed = 1,
                              sample.size = 0.9 * min(sample_sums(ps)))
# compute phylogenetic diversity
pd <-
  metagMisc::phyloseq_phylo_div(ps_rarefied) %>%
  tibble::rownames_to_column("sample_name")

# bind metadata
pd_meta <-
  dplyr::left_join(
    data.frame(sample_data(ps_rarefied)),
    pd,
    by = "sample_name")

#plot data
ggplot(pd_meta) +
  geom_point(aes(x = species, y = PD, color = season)) +
  facet_wrap(~rhizosphere)
ggsave(file = "output/phylogenetic_diversity.pdf")
