###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: visualize trees with metabarcoder for significant contrasts:
# Spartina maritima: rhizosphere vs bulk-soil - check groups for S. maritima
# species: all against Spartina maritima - check for groups related to low marsh
#PROJECT: spartina-metarizo
###.............................................................................
library(metacoder)
library(phyloseq)
library(tidyverse)

# read phyloseq
ps <-
  readRDS("data/intermediate/ps_t_noRep.rds") %>%
  phyloseq::subset_samples(species == "Spartina_maritima")


# convert phyloseq to taxmap
ps_taxmap <- metacoder::parse_phyloseq(ps)

# compute otu_table with proportions. Alternatively, data can be rarefied/normalized.
ps_taxmap$data$tax_data <-
  calc_obs_props(ps_taxmap, "otu_table")

# compute abundance per taxa
ps_taxmap$data$tax_abund <-
  calc_taxon_abund(ps_taxmap,
                   "tax_data", # otu_table with proportions
                   cols = ps_taxmap$data$sample_data$sample_id) # per sample

# compute ocurrence of each taxon per sample_data category
ps_taxmap$data$tax_occ <-
  metacoder::calc_n_samples(ps_taxmap,
                            "tax_abund",
                            groups = ps_taxmap$data$sample_data$rhizosphere,
                            cols = ps_taxmap$data$sample_data$sample_id)

# compare grups
ps_taxmap$data$diff_table <- metacoder::compare_groups(ps_taxmap,
                                      dataset = "tax_abund", # proportions
                                      cols = ps_taxmap$data$sample_data$sample_id, # What columns of sample data to use
                                      groups = ps_taxmap$data$sample_data$rhizosphere) # What category each sample is assigned to
# plot comparissons
set.seed(1)
heat_tree_matrix(ps_taxmap,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree_smrizo.pdf") # Saves the plot as a pdf file
