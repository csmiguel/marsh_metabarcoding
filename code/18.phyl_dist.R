###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: turnover phylogenetic diversity
#PROJECT: spartina-metarizo
###.............................................................................
library(reshape2)
library(picante)
library(phyloseq)
library(tidyverse)


nperm <- 2 # number of permutations

source("code/functions/ps_filter_prevalence.r") # function to filter phyloseq
source("code/functions/bmntd.r") # function to filter phyloseq

ps <- readRDS("data/intermediate/ps_t_noRep.rds") # load phyloseq

##0bserved values for mean nearest taxon distance
otu_obs <- as.data.frame(otu_table(ps)) # observed otu table
ctree_obs <- cophenetic(phy_tree(ps)) # cophenetic distance
bmntd_obs <- # observed B mean nearest taxon distance
  picante::comdistnt(comm = otu_obs,
                   dis = ctree_obs,
                   abundance.weighted = T,
                   exclude.conspecifics = T)
bmntd_obs <- # distance object to data frame
  melt(as.matrix(bmntd_obs), varnames = c("sample1", "sample2"))
### 1 ### phylogenetic filtering in rhizosphere
# generate data frame with contrasts to test; one per row
#  this set of contrast will test phylogenetic turnover between bulk-> rhizospheric soil
contrasts_rhiz <-
  as(sample_data(ps), "data.frame") %>%
  mutate(selectioncol = paste0(species, season)) %>%
  filter(duplicated(selectioncol)) %>%
  select(species, season) %>%
  arrange(species, season)
# run permutations
rhiz_perm <-
  contrasts_rhiz %>%
  apply(1, function(x) { # per line
    sp <<- x[1]
    sea <<- x[2]
    # subset phyloseq according to contrast
    ps1 <-
      subset_samples(ps, species == sp & season == sea) %>%
      ps_filter_prevalence(prevalence_threshold = 1) %>%
      transform_sample_counts(function(x) x / sum(x))
    phydistperm <-
      permutate_bmntd(physeq = ps1, nperm = nperm)
    phydistperm %>%
      mutate(contrast = paste(sp, sea, sep = "_"))
  })
# df with summary statistics
rhiz_perm_summ <-
  rhiz_perm %>%
    lapply(function(x){
      data.frame(contrast = unique(x$contrast),
                 sample1 = unique(x$row),
                 sample2 = unique(x$col),
                 sdev = sd(x$value),
                 av = mean(x$value))
    }) %>%
    do.call(what = "rbind")

### 2 ### phylogenetic filtering in between seasons
# generate data frame with contrasts to test; one per row
# this set of contrast will test phylogenetic turnover between subsecuent seasons for bulk soil
contrasts_season <-
  as(sample_data(ps), "data.frame") %>%
  filter(rhizosphere == "bulk_soil") %>%
  select(sample_name, species, season) %>%
  dcast(species ~ season, value.var = "season") %>%
  mutate(season1 = paste(spring, summer),
         season2 = paste(summer, autumn),
         season3 = paste(autumn, winter)) %>%
  melt(id.vars = "species", measure.vars =  c("season1", "season2", "season3")) %>%
  filter(!grepl("NA", value)) %>%
  tidyr::separate(col = value, c("season1", "season2"), sep = " ") %>%
  select(-variable)

# run permutations
season_perm <-
  contrasts_season %>%
  apply(1, function(x) { # per line
    sp <<- x[1]
    sea1 <<- x[2]
    sea2 <<- x[3]
    # subset phyloseq according to contrast
    ps1 <-
      subset_samples(ps, species == sp & rhizosphere == "bulk_soil") %>%
      subset_samples(season == sea1 | season == sea2) %>%
      ps_filter_prevalence(prevalence_threshold = 1) %>%
      transform_sample_counts(function(x) x / sum(x))
    phydistperm <-
      permutate_bmntd(physeq = ps1, nperm = nperm)
    phydistperm %>%
      mutate(contrast = paste(sp, sea1, sea2, sep = "_"))
  })
# df with summary statistics
season_perm_summ <-
  season_perm %>%
  lapply(function(x){
    data.frame(contrast = unique(x$contrast),
               sample1 = unique(x$row),
               sample2 = unique(x$col),
               sdev = sd(x$value),
               av = mean(x$value))
  }) %>%
  do.call(what = "rbind")

# merge tables
rbind(rhiz_perm_summ, season_perm_summ) %>%
  left_join(x = bmntd_obs,
            by = c("sample1", "sample2")) %>%
  filter(!is.na(sdev))
