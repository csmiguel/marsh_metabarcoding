###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: infere taxa with differential abundances by selected variables
#PROJECT: spartina-metarizo
###...........................................................................

library(tidyverse)
library(DESeq2)
library(phyloseq)

# load functions
source("code/functions/gm_means.r") # geometric means
source("code/functions/deseqp.r") # table from deseq2 result object
source("code/functions/deseq2_tests.r") # function to carry tests

## 0. agglomarate phyloseq according to taxonomic rank
# set rank to perform analysis
ranki <- "family"

# read phyloseq object and agglomarate counts according to rank
ps <-
  readRDS("data/intermediate/ps_t_noRep.rds") %>%
  phyloseq::tax_glom(ranki)

## test 1: differentianl abundaces between species plot
# all crossed factors
# unique combinations for species
combination_species <-
  unique(sample_data(ps)$species) %>%
    combn(2, simplify = F)
# tests
tests_species <-
    lapply(combination_species, function(x) c("species", x))
# deseq test
ds_species <- ds_test(ps = ps, tests_species)

## test 2: differentianl abundaces between rhizosphere vs bulk soil
#   tests
tests_rhizosphere <- list(c("rhizosphere", "rhizosphere", "bulk_soil"))
#   deseq test
ds_rhizosphere <- ds_test(ps = ps, tests_rhizosphere)

## test 3: differentianl abundaces between rhizosphere vs bulk soil for each species
#   split phyloseq per species
sp <- unique(sample_data(ps)$species)
#   deseq test
ds_rhizosphere_sp <-
  sp %>%
    lapply(function(x) {
      temp_x <<- x # for some reason functions in phyloseq fail to be fed with
                  # objects from the function env
      #subset ps per species
        phyloseq::subset_samples(ps, species == temp_x) %>%
         #rm ASVs with 0 abundances
        metagMisc::phyloseq_filter_taxa_rel_abund(0) %>%
        ds_test(tests_rhizosphere)
    })

## test 4: differentianl abundaces across seasons
#   unique combinations for species
combination_seasons <- combn(levels(sample_data(ps)$season), 2, simplify = F)

#   tests
tests_seasons <-
    lapply(combination_seasons, function(x) c("season", x))
#   deseq test
ds_seasons <- ds_test(ps = ps, tests_seasons)

## part2: generate taxa dataframe with agglomerated taxa by family
# melt phyloseq
psmelted <- phyloseq::psmelt(ps)

## 1. create data frame with ASV names and taxa names
# formula to cast melted phyloseq
fcast <-
  rank_names(ps)[1:match(ranki, rank_names(ps))] %>%
      paste(collapse = " + ") %>%
      paste("OTU +", ., "~ sample_name")

#tax data frame
tax_df <-
  reshape2::dcast(psmelted, fcast,
    value.var = "Abundance") %>%
    dplyr::select("OTU", rank_names(ps)[1:match(ranki, rank_names(ps))])

#save objects
saveRDS(
  list(test1 = ds_species,
        test2 = ds_rhizosphere,
        test3 = ds_rhizosphere_sp,
        test4 = ds_seasons),
      "data/intermediate/deseq2_results.rds")
saveRDS(tax_df, "data/intermediate/tax_df.rds")
