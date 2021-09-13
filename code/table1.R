#produce table 1

library(phyloseq)
library(tidyverse)
library(flextable)

# extract metadata from phyloseq object
meta <-
  readRDS("data/intermediate/ps_marsh_filt.rds") %>%
  phyloseq::sample_data() %>% # select table with metadata
  as("data.frame") %>%
  dplyr::filter(rep < 2) %>% # remove replica
  dplyr::select(sample_name,
                species,
                season,
                rhizosphere) %>% # select variables of interest
  dplyr::arrange(species, season, rhizosphere) %>%
  dplyr::rename("ID" = sample_name,
                "Species plot" = species,
                "Season" = season,
                "Soil" = rhizosphere) %>%
  as_tibble()

# sample sequencing depth
seq_depth <-
  read.table("output/sequencing_report.txt", header = T) %>%
  dplyr::filter(file != "file") %>%
  dplyr::mutate(sample = stringr::str_extract(file, "S-[0-9]++"),
                depth = as.numeric(stringr::str_remove(num_seqs, ","))) %>%
  dplyr::select(sample, depth)

# format table
ft <-
  flextable(meta) %>%
  set_caption(caption = "Samples") %>%
  italic(j = 2) %>%
  bold(part = "header", bold = TRUE)
ft
