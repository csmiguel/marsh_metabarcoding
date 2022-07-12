###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: phylum plots
#PROJECT: marsh_metabarcoding
###.............................................................................
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)

# load phyloseq
ps <-
  readRDS("data/intermediate/ps_t_noRep.rds.rds") %>%
  phyloseq::transform_sample_counts(function(x) { # transform to proportions
     x / sum(x) }) %>%
  phyloseq::tax_glom(taxrank = "phylum", NArm = FALSE) # agglomerate

#get top phyla with more abundances
names_top_phyla <-
  colSums(otu_table(ps)) %>%
  sort(decreasing = T) %>%
  names() %>%
  .[1:10]

#tax_table
tax_table_phylum <-
  tax_table(ps)@.Data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("asv")

#vector with most abundant phyla
vector_names_legend <-
  c(tax_table_phylum$phylum[match(names_top_phyla, tax_table_phylum$asv)],
    "other")

#get tidy counts of abundance per phylum
tidy_counts_phylum <-
  phyloseq::psmelt(ps) %>%
  as_tibble() %>%
  group_by(sample_name, OTU, phylum, sample_species, rhizosphere, season) %>%
  rename(species = sample_species) %>%
  summarise(sum_phylum = sum(Abundance)) %>%
  #rename to "other" all phyla not with abundances
  mutate(to_phyla = ifelse(OTU %in% names_top_phyla,
                           phylum,
                           "other")) %>%
  group_by(sample_name, to_phyla, species, rhizosphere, season) %>%
  summarise(sum_phylum = sum(sum_phylum)) %>%
  mutate(Phylum = factor(to_phyla,
                         levels = vector_names_legend),
        species = gsub("_", "\n", species),
        rhizosphere = gsub("_soil", " ", rhizosphere))

#barplot
p1 <-
  ggplot(tidy_counts_phylum) +
  geom_bar(aes(x = season, y = sum_phylum, fill = Phylum),
           stat = "identity",
           position = "stack",
           color = "black",
           size = 0.1) +
  facet_grid(rhizosphere ~ species, scales = "fixed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "italic")) +
  scale_fill_brewer(palette = "Paired") +
  ylab("Relative abundance")

#save plot
prop_p1 <- 9 / 7
ggsave(filename = "output/phylum-barplot.pdf", p1,
        width = prop_p1 * 6, height = 6)
