###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: beta diversity: MDS plots from unifrac distances
#PROJECT: spartina-metarizo
###.............................................................................
library(tidyverse)
library(phyloseq)
library(ggrepel)
library(cowplot)

ps <-
  readRDS("data/intermediate/ps_t.rds") %>%
  phyloseq::subset_samples(rep < 2)

# transformation of the abundances to natural logarithmic scale
ps_log <- phyloseq::transform_sample_counts(ps, function(x) log(1 + x))

#create grouping factor to connect rhizosphere/bulk soil samples in plot
sample_data(ps_log) <-
  sample_data(ps_log) %>%
    as("data.frame") %>%
    dplyr::mutate(groupr = paste0(season, species)) %>%
    sample_data()

# MDS with weighted UNIFRAC distances
unif_meta <-
  phyloseq::ordinate(ps_log,
    method = "MDS",
    distance = "wunifrac")

#get eigenvalues
evals <- unif_meta$values$Eigenvalues

#plot MDS
# PC1 vs PC2
p1 <-
  phyloseq::plot_ordination(ps_log,
                            unif_meta,
                            axes = c(1, 2),
                            shape = "species",
                            color = "season") +
    labs("species", "season") +
    geom_line(aes(group = groupr), color = "grey") +
    geom_point(size = 2) +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    ggrepel::geom_text_repel(
      aes(label = sample_data(ps_log)$rhizosphere %>% substring(1, 1)),
      size = 2,
      min.segment.length = 1,
      point.padding = 0.5) +
    scale_colour_manual(
      values = c("darkgreen", "orange", "brown", "darkblue")) +
    theme_cowplot()

# PC3 vs PC2
p2 <-
  phyloseq::plot_ordination(ps_log,
                            unif_meta,
                            axes = c(3, 2),
                            shape = "species",
                            color = "season") +
  labs("species", "season") +
  geom_line(aes(group = groupr), color = "grey") +
  geom_point(size = 2) +
  coord_fixed(sqrt(evals[2] / evals[3])) +
  ggrepel::geom_text_repel(
    aes(label = sample_data(ps_log)$rhizosphere %>% substring(1, 1)),
    size = 2,
    min.segment.length = 1,
    point.padding = 0.5) +
  scale_colour_manual(
    values = c("darkgreen", "orange", "brown", "darkblue")) +
  theme_cowplot()

# plot dimensions
pcs <-
  data.frame(value = 100 * unif_meta$values$Relative_eig,
             Dimension = seq_along(unif_meta$values$Relative_eig)) %>%
  head(n = 8)

p_dims <-
  ggplot2::ggplot(aes(x = Dimension, y = value),
                data = pcs) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = pcs$Dimension, breaks = pcs$Dimension) +
  xlab("Dimension") +
  ylab("Variation %") +
  theme_cowplot()

# composite plot
top_row <-
  plot_grid(p1)
bottom_row <-
  plot_grid(p2 + theme(legend.position = "none"),
            p_dims)
pall <-
  plot_grid(top_row,
          bottom_row,
          ncol = 1)
# facet plot
# 1vs2
p_facet1vs2 <-
  phyloseq::plot_ordination(ps_log,
                            unif_meta,
                            axes = c(1, 2),
                            color = "season") +
  labs("species", "season") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_line(aes(group = groupr), color = "black", lwd = 0.5) +
  geom_point(size = 2) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  ggrepel::geom_text_repel(
    aes(label = sample_data(ps_log)$rhizosphere %>% substring(1, 1)),
    size = 3,
    min.segment.length = 1,
    point.padding = 0.5) +
  scale_colour_manual(
    values = c("darkgreen", "orange", "brown", "darkblue")) +
  facet_wrap(~species) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey80",
                                        linetype = "dotted",
                                        size = 0.3))
# 3vs4
p_facet3vs4 <-
  phyloseq::plot_ordination(ps_log,
                            unif_meta,
                            axes = c(3, 4),
                            color = "season") +
  labs("species", "season") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_line(aes(group = groupr), color = "black", lwd = 0.5) +
  geom_point(size = 2) +
  coord_fixed(sqrt(evals[4] / evals[3])) +
  ggrepel::geom_text_repel(
    aes(label = sample_data(ps_log)$rhizosphere %>% substring(1, 1)),
    size = 3,
    min.segment.length = 1,
    point.padding = 0.5) +
  scale_colour_manual(
    values = c("darkgreen", "orange", "brown", "darkblue")) +
  facet_wrap(~species) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey80",
                                        linetype = "dotted",
                                        size = 0.3))

#save plots
ggsave("output/MDS_unifrac.pdf",
       pall,
       width = 9,
       height = 7,
       units = "in")

ggsave("output/MDS_pfacet1vs2.pdf",
       p_facet1vs2,
       width = 9,
       height = 5,
       units = "in")

ggsave("output/MDS_pfacet3vs4.pdf",
      p_facet3vs4,
      width = 9,
      height = 5,
      units = "in")

#save transformed counts for weighted unifrac distace calculations
saveRDS(ps_log, "data/intermediate/ps_log.rds")
