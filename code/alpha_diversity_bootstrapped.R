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
library(reshape2)

#read filtered phyloseq
ps <- readRDS("data/intermediate/ps_filt.rds")
#load functions
source("code/functions/alpha_diversity_bt.r")
source("code/functions/alpha_diversity.r")

alpha_outliers <- alpha_diversity(ps = ps, "output/alpha_outliers_vs_readcount.pdf")
#alfter visual inspection of the output/alpha_outliers.pdf I remove ones sample:
ps_pruned <- phyloseq::prune_samples(sample_names(ps)[!sample_names(ps) %in% "S-16"], ps)

# estimate alpha diversity over bootstrapped OTU matrices
adiv_bt <-
  alpha_diversity_bt(ps = ps_pruned, bt = 100)
adiv_bt2 <-
  adiv_bt %>%
  dplyr::filter(rep == 1) %>%
  dplyr::select(-dataset, -co2, -temp, -salt, -rep) %>%
  reshape2::melt(id.vars = c("sample_name", "species", "season", "rhizosphere"),
                 value.name = "Shannon") %>%
  dplyr::mutate(season = as.factor(season)) %>%
  as_tibble()

#plot bootstrapped alpha diversity
p1 <-
  adiv_bt2 %>%
  ggplot() +
  geom_boxplot(aes(x = season,
                 y = Shannon, color = rhizosphere),
               outlier.shape = NA,
               position = "identity",
               width = 0.1) +
  geom_line(data = adiv_bt2 %>%
              group_by(species, season, rhizosphere) %>%
              summarise(mean_shannon = mean(Shannon)),
            aes(x = as.numeric(season),
                y = mean_shannon, linetype = rhizosphere)) +
  facet_grid(~species) +
  #scale_x_discrete(limits = levels(as.factor(data_glm$season))) +
  theme_classic() +
  scale_color_manual(values = c("black", "black")) +
  theme(
    panel.grid.major =
      element_line(colour = "grey50",
                   size = 0.1),
    axis.text.x = element_text(angle = 90),
    strip.background = element_blank()) +
  guides(colour = F)
#write plot
ggsave(filename = "output/alpha-diversity_bt.pdf", p1, width = 12, height = 3)
