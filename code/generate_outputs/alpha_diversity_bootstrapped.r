# plot alpha diversity bootstrapped
library(ggplot2)
library(dplyr)

data <-
  readRDS("data/intermediate/alpha_div_bt.rds")

#plot bootstrapped alpha diversity
p1 <-
  ggplot(data) +
  geom_boxplot(aes(x = season,
                 y = Shannon, color = rhizosphere),
               outlier.shape = NA,
               position = "identity",
               width = 0.1) +
  geom_line(data = data %>%
              group_by(species, season, rhizosphere) %>%
              summarise(mean_shannon = mean(Shannon)),
            aes(x = as.numeric(season),
                y = mean_shannon, linetype = rhizosphere)) +
  facet_grid(~species) +
  theme_classic() +
  scale_color_manual(values = c("black", "black")) +
  ylab("Shannon H'") +
  theme(
    panel.grid.major =
      element_line(colour = "grey50",
                   size = 0.1),
    axis.text.x = element_text(angle = 90),
    strip.background = element_blank()) +
  guides(colour = F)
#write plot
ggsave(filename = "output/alpha-diversity_bt.pdf", p1, width = 12, height = 3)
