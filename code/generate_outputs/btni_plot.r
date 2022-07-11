library(tidyverse)
library(reshape2)

bnti <-
  readRDS("data/intermediate/bnti.rds") %>%
  tidyr::separate(contrast, c("genus", "sp", "var1", "var2"), sep = "_") %>%
  mutate(species = paste(genus, sp, sep = "_")) %>%
  select(-genus, -sp) %>%
  as_tibble()

# filtro rhizosfera
rhizo <-
  filter(bnti, var2 == "bulk2rhizo") %>%
  select(-sample1, -sample2, -var2) %>%
  rename(season = var1) %>%
  mutate(season = factor(season, levels = c("spring", "summer", "autumn", "winter")),
         significant = as.factor(abs(bnti) > 2))

p_rhizo <-
  ggplot(rhizo, aes(x = season, y = observed -av)) +
  geom_point(aes(color = significant), size = abs(rhizo$bnti)) +
  geom_line(aes(as.numeric(as.factor(season))), linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = observed -av - sdev, ymax = observed -av + sdev),
                width = 0.5) + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  facet_wrap(~species, nrow = 1, scales = "free_x") +
  scale_color_manual(values = c("black", "red")) +
  ylab("bMNTDobs - bMNTDnull") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background = element_blank())

ggsave("output/bnti_rhizo.pdf", width = 11, height = 4)

# filtro season
sesonturn <-
  filter(bnti, var2 != "bulk2rhizo") %>%
  select(-sample1, -sample2) %>%
  mutate(season = factor(paste(var1, var2, sep = "_"),
                         levels = c("spring_summer", "summer_autumn", "autumn_winter")),
         significant = as.factor(abs(bnti) > 2))

# plot
p_season <-
  ggplot(sesonturn, aes(x = season, y = observed - av)) +
  geom_point(aes(color = significant), size = abs(sesonturn$bnti)) +
  geom_line(aes(as.numeric(as.factor(season))), linetype = "dashed", alpha = 0.5) +
  geom_errorbar(aes(ymin = observed -av - sdev, ymax = observed -av + sdev),
                width = 0.5) + 
  geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
  facet_wrap(~species, nrow = 1, scales = "free_x") +
  ylab("bMNTDobs - bMNTDnull") +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background = element_blank())
ggsave("output/bnti_season.pdf", width = 9, height = 4)
