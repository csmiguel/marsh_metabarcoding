###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: beta diversity: dendrogram from wnifrac
#PROJECT: spartina-metarizo
###.............................................................................
library(tidyverse)
library(phyloseq)

ps <- readRDS("data/intermediate/ps_log.rds")

# weighted wunifrac distances
wuni <- phyloseq::distance(ps, method = "wunifrac")

# hierarchical clustering
hc <- hclust(wuni, method = "complete")
hc_phylo <- ape::as.phylo(hc)

#data frame with metadata for mapping to the tree
meta <-
  sample_data(ps) %>% as("data.frame") %>%
    select(sample_name, rhizosphere, species, season)

#plot
p <-
  ggtree(hc_phylo) %<+% meta +
  geom_tippoint(
        mapping = aes(x = x + 0.0005, shape = species),
        size = 2) +
        geom_tiplab(aes(geom = "sample_name"), offset = 0.001, size = 2) +
        geom_treescale(label = "weighted unifrac distance")

# data frame for mapping heatmap
meta2 <- meta[, 4, drop = F]

# add heatmap
p2 <-
  gheatmap(p, meta2, offset = 0.003, width = 0.05, hjust = 0) +
  scale_fill_manual(
    values = c("darkgreen", "orange", "brown", "darkblue")) +
  scale_shape_manual(values = c(0, 1, 2, 7, 9))

#save plot
ggsave("output/dendrogram_wunifrac.pdf",
       p2,
       width = 9,
       height = 7.7,
       units = "in")
