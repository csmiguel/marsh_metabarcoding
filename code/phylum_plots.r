###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: phylum plots
#PROJECT: spartina-metarizo
###.............................................................................
library(ggplot2)
library(dplyr)
library(phyloseq)

ps <- readRDS("data/intermediate/ps_filt.rds")
#rarefy
ps_filt <-
  phyloseq::prune_samples(
    sample_names(ps)[!sample_names(ps) %in% "S-16"], ps) %>%
  phyloseq::subset_samples(rep == 1) %>%
  phyloseq::transform_sample_counts(function(x) { x / sum(x) })

#agglomerate by phylum
ps_phylum <-
  tax_glom(ps_filt,
           taxrank = "phylum",
           NArm = FALSE)

#get top phyla with more abundances
names_top_phyla <-
  colSums(otu_table(ps_phylum)) %>%
  sort(decreasing = T) %>%
  names() %>%
  .[1:10]
#tax_table
tax_table_phylum <-
  tax_table(ps_phylum)@.Data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("asv")

#vector with most abundant phyla
vector_names_legend <-
  c(tax_table_phylum$phylum[match(names_top_phyla, tax_table_phylum$asv)], "other")
#get tidy counts of abundance per phylum
tidy_counts_phylum <-
  phyloseq::psmelt(ps_phylum) %>%
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
                         levels = vector_names_legend))

#barplot
p1 <-
  ggplot(tidy_counts_phylum) +
  geom_bar(aes(x = season, y = sum_phylum, fill = Phylum),
           stat = "identity", position = "stack", color = "black") +
  facet_grid(rhizosphere~species, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0),
        strip.background = element_blank()) +
  ylab("Proportion of counts")

#save plot
ggsave(filename = "output/phylum-barplot.pdf", p1, width = 9, height = 7)
