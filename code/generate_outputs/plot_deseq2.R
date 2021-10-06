###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# oct 2021
###.............................................................................
#GOAL: plot deseq results per factor
#PROJECT: marsh_metabarcoding
###.............................................................................

library(tidyverse)
library(phyloseq)
library(pheatmap)

# read tax table
tax <-
  readRDS("data/intermediate/tax_df.rds") %>%
  mutate(phyfam = paste(phylum, family, sep = ":")) %>%
  dplyr::select("OTU", "phyfam")
# > head(tax)
#     OTU                                phyfam
# 1  ASV1        Campilobacterota:Sulfurovaceae
# 2  ASV9                Firmicutes:Bacillaceae
# 3  ASV4 Proteobacteria:Pseudoalteromonadaceae
# 4 ASV18   Desulfobacterota:Desulfosarcinaceae
# 5 ASV60           Proteobacteria:Vibrionaceae
# 6 ASV62  Actinobacteriota:Geodermatophilaceae

# read deseq objects
ds <- readRDS("data/intermediate/deseq2_results.rds")

# season
plot_deseq2(ds_result = ds$test4, out_file = "output/deseq2_season.pdf")

# species
plot_deseq2(ds_result = ds$test1, out_file = "output/deseq2_species.pdf")

# rhizosphere
plot_deseq2(ds_result = c(ds$test2, ds$test3),
            out_file = "output/deseq2_rhizosphere.pdf", use_id = T)
