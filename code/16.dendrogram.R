###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# August 2021
###.............................................................................
#GOAL: beta diversity: dendrogram from wnifrac
#PROJECT: marsh_metabarcoding
###.............................................................................
library(dplyr)
library(phyloseq)
library(ggtree)
library(ggplot2)

ps <-
  readRDS("data/intermediate/ps_t_noRep.rds") %>%
  # transformation of the abundances to natural logarithmic scale
  phyloseq::transform_sample_counts(function(x) log(1 + x))

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
#be aware of the following issue: https://github.com/YuLab-SMU/ggtree/issues/399
# and how to solve it:
#require(devtools)
#install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")
p <-
  ggtree(hc_phylo) %<+% meta +
  geom_tippoint(
        mapping = aes(x = x + 0.0005, shape = species, color = season),
        size = 2) +
        geom_tiplab(aes(geom = "sample_name"), offset = 0.001, size = 2) +
        geom_treescale(label = "weighted unifrac distance") +
        scale_color_manual(
          values = c("darkgreen", "orange", "brown", "darkblue")) +
        scale_shape_manual(values = c(0, 1, 2, 7, 9))

#save plot
ggsave("output/dendrogram_wunifrac.pdf",
       p,
       width = 9,
       height = 7.7,
       units = "in")
