###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: compute alpha diversity
#PROJECT: marsh_metabarcoding
###.............................................................................
library(phyloseq)
library(ggplot2)
library(dplyr)

#read filtered phyloseq
ps <- readRDS("data/intermediate/ps_t_noRep.rds")

#load functions
source("code/functions/ps2df.r")
source("code/functions/alpha_diversity.r")
source("code/functions/alpha_diversity_bt.r")

#estimate alpha diversity
alpha_div <-
  alpha_diversity(ps = ps, "output/alpha_vs_readcount.pdf") %>%
  tibble::rownames_to_column("sample_name")

# estimate Shannon diversity over bootstrapped OTU matrices
alpha_div_bt <-
  alpha_diversity_bt(ps = ps, bt = 100) %>%
  dplyr::select(-rep) %>%
  reshape2::melt(id.vars = c("sample_name", "species", "season", "rhizosphere"),
                 value.name = "Shannon") %>%
  as_tibble()

# save objects
saveRDS(alpha_div, "data/intermediate/alpha_div.rds")
saveRDS(alpha_div_bt, "data/intermediate/alpha_div_bt.rds")

# write alpha diversity indices to output

write.csv(alpha_div, "output/alpha_diversity.csv")
