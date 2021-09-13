###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: model change of alpha diversity across conditions
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(ggplot2)
library(lme4)
library(dplyr)

#read filtered phyloseq
ps <- readRDS("data/intermediate/ps_filt.rds")
#load functions
source("code/functions/ps2df.r")
source("code/functions/alpha_diversity.r")

#alfter visual inspection of the output/alpha_outliers.pdf I remove ones sample:
ps_pruned <-
  phyloseq::prune_samples(sample_names(ps)[!sample_names(ps) %in% "S-16"], ps)

# alpha diversity
div <-
  ps_pruned %>%
    phyloseq::estimate_richness(measures = "Shannon") %>%
    tibble::rownames_to_column("sample_name")

#tidy data for models
data_glm <-
  ps2df(ps_pruned)$sample_data %>%
    dplyr::left_join(div, by = "sample_name") %>%
    dplyr::select(-co2, -temp, -salt) %>% as_tibble()

#check
#https://ourcodingclub.github.io/tutorials/mixed-models/#six
# Is diversity dependent on the rhizosphere?
#div ~rhizosphere: the
m1 <- lme4::lmer(Shannon ~ species:rhizosphere + (1|season), data = data_glm)
summary(m1)
ranef(m1)
sjPlot::plot_model(m1,
                   #axis.labels=c("Urchin", "Depth", "Fish"),
                   show.values = TRUE, show.p = TRUE,
                   title = "Effect of rhizosphere on diversity")
# Is diversity dependent on the season? div ~season
m2 <- lme4::lmer(Shannon ~ season + (1|species:rhizosphere), data = data_glm)
summary(m2)
sjPlot::plot_model(m2,
                   #axis.labels=c("Urchin", "Depth", "Fish"),
                   show.values = TRUE,
                   show.p = TRUE,
                   sort.est = TRUE, type = "est",
                   show.intercept = T,
                   title = "Effect of season on diversity")
sjPlot:: tab_model(m2)
# Is diversity dependent on the species?
m3 <- lme4::lmer(Shannon ~ species:rhizosphere + (1|season), data = data_glm)
summary(m3)
ranef(m3)
sjPlot::plot_model(m3,
                   #axis.labels=c("Urchin", "Depth", "Fish"),
                   show.values = TRUE, show.p = TRUE,
                   title = "Effect of rhizosphere on diversity")


#prepare data for plotting
h <-
  get_model_data(model = m2, type = "est",
               show.intercept = T,
               ci.lvl = 0.95) %>%
  as.data.frame() %>%
  dplyr::select(1:3) %>%
  dplyr::mutate(season = levels(as.factor(data_glm$season))) %>%
  dplyr::select(season, estimate, std.error)
h[1, 2] <- 0

ggplot(h, aes(x = season, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = .1) +
  ggtitle("Effect of season on alpha diversity")
