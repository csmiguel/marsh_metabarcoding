###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: compute alpha diversity bootstrapped
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(ggplot2)
library(dplyr)

#read filtered phyloseq
ps <- readRDS("data/intermediate/ps_filt.rds")
#load functions
source("code/functions/ps2df.r")
source("code/functions/alpha_diversity.r")

alpha_outliers <-
  alpha_diversity(ps = ps, "output/alpha_outliers_vs_readcount.pdf")
#alfter visual inspection of the output/alpha_outliers.pdf I remove ones sample:
ps_pruned <-
  phyloseq::prune_samples(sample_names(ps)[!sample_names(ps) %in% "S-16"], ps)

# and recalculate alpha diversity
alpha_div <-
  alpha_diversity(ps_pruned, "output/alpha-diversity_vs_readcount.pdf") %>%
  tibble::rownames_to_column("sample_name")

# change of alpha diversity with season
p1 <-
  data_glm %>% as_tibble() %>%
  mutate(season = as.numeric(as.factor(season))) %>%
  filter(rep == 1) %>%
  ggplot() +
  geom_line(aes(x = season,
                y = Shannon,
                linetype = rhizosphere)) +
  facet_grid(~species) +
  scale_x_discrete(limits = levels(as.factor(data_glm$season))) +
  theme_classic() +
  theme(
    panel.grid.major =
      element_line(colour = "grey50",
                   size = 0.1),
    axis.text.x = element_text(angle = 90),
    strip.background = element_blank())
#write plot
ggsave(filename = "output/alpha-diversity.pdf", p1, width = 12, height = 3)
