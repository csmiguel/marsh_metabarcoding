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
library(lme4)
library(tidyverse)
library(car)
library(sjPlot)

# read diversity
div <- readRDS("data/intermediate/alpha_div.rds")

# read sample metadata
meta <-
  readRDS("data/intermediate/ps_t_noRep.rds") %>%
  phyloseq::sample_data() %>%
  as("data.frame") %>%
  dplyr::select(-rep)

#tidy data for models
data <-
    dplyr::left_join(div, meta, by = "sample_name") %>%
    as_tibble()
# check https://strengejacke.github.io/sjPlot/articles/tab_model_estimates.html
# Is diversity dependent on the season? div ~season
m11 <- lme4::lmer(Shannon ~ season + (1|species), data = data)
m10 <- lmer(Shannon ~ (1|species), data = data)
summary(m11)
t1 <- sjPlot:: tab_model(m11)
anova(m11, m10)

# Is diversity dependent on the rhizosphere for each species?
m21 <-
  lme4::lmer(Shannon ~ species * rhizosphere + (1|season), data = data)
summary(m21)
m20 <- lme4::lmer(Shannon ~ species + rhizosphere + (1|season), data = data)
t2 <- sjPlot:: tab_model(m21)
anova(m21, m20)

#output table
tab_model(m11, m21,
          dv.labels = c("Shannon ~ season", "Shannon ~ species * rhizosphere"),
        file = "output/alpha_models.doc")

# do different species have diferent variance in diversity across seasons?
# levene test for inequlity of variance for each four data points of soil within
# each of the species
levtest <- car::leveneTest(Shannon ~ species:rhizosphere, data = data)
# create dataset with sd and means per group
d1 <-
  data %>%
    dplyr::group_by(species, rhizosphere) %>% # Specify group indicator
    dplyr::summarise_at(vars(Shannon),        # Specify column
                   list(sd = sd, mean = mean)) %>%
    dplyr::mutate(sd = round(sd, 2))
# save results to file
save(levtest, d1, "data/intermediate/levene.Rdata")
