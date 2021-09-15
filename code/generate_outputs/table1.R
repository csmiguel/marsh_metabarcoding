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
              "No. reads" = as.numeric(stringr::str_remove(num_seqs, ","))) %>%
  dplyr::select(sample, "No. reads") %>%
  dplyr::rename(ID = sample)

# no of ASVs per sample
asvs <-
  readRDS("data/intermediate/ps_filt.rds") %>%
  otu_table() %>%
  apply(1, function(x) sum(x > 0)) %>%
  data.frame(ID = names(.), ASVs = .)

# join data
data <-
  dplyr::left_join(meta, seq_depth, by = "ID") %>%
    dplyr::left_join(asvs, by = "ID") %>%
    dplyr::mutate("Species plot" = gsub("_", " ", `Species plot`),
                  Soil = gsub("_soil", "", Soil))

# format table
vector_borders <- cumsum(data$`Species plot` %>% table %>% as.numeric())
ft <-
  flextable(data) %>%
  set_caption(caption = "Table 1. Sample characteristics, sequencing depth
 (No. reads), number of ASVs and SRA accession.") %>%
  italic(j = 2) %>%
  fontsize(part = "all", size = 10) %>%
  merge_v(j = c("Species plot", "Season")) %>%
  set_table_properties(width = .5, layout = "autofit") %>%
  bold(part = "header", bold = TRUE) %>%
  hline(i = vector_borders,
        part = "body") %>%
  fix_border_issues()

# save to docx
save_as_docx(ft, path = "output/table1.docx")
