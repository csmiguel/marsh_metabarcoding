###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: beta diversity. ADONIS permutation test
#PROJECT: spartina-metarizo
###.............................................................................
library(vegan)
library(phyloseq)
library(dplyr)
library(flextable)

# variable with permutations
nperm <- 9999
set.seed(123)

# read phyloseq
ps <- readRDS("data/intermediate/ps_t_noRep.rds")  %>%
# transformation of the abundances to natural logarithmic scale
phyloseq::transform_sample_counts(function(x) log(1 + x))

#data frame with metadata
metadf <- data.frame(sample_data(ps))

# weighted unifrac distances
wudist <- phyloseq::UniFrac(ps,
                    weighted = TRUE,
                    normalized = TRUE,
                    parallel = FALSE,
                    fast = TRUE)

# test 1: effect of species
t1 <-
  vegan::adonis(
    wudist ~ species,
    data = metadf,
    permutations = nperm)

# test 2: effect of rhizosphere
# across all samples and seasons
t2 <-
  vegan::adonis(
    wudist ~ rhizosphere,
    data = metadf,
    permutations = nperm)
# is marginal significant, but as we saw in the MDS the variation of the
#   community has different directions depending mainly on the species, and
#   partly on the season.

# test 3: effect of rhizosphere per species
t3 <-
  vegan::adonis(
    wudist ~ rhizosphere,
    data = metadf,
    permutations = nperm,
    strata = metadf$species)

# test 4: change with season
t4 <-
  vegan::adonis(
    wudist ~ season,
    permutations = nperm,
    data = metadf)

# test 5: change with season for every species
t5 <-
  vegan::adonis(
    wudist ~ season,
    data = metadf,
    permutations = nperm,
    strata = metadf$species)

# report results from adonis
sink(file = "output/adonis.txt")
cat("Results from vegan::adonis analysis:
  test1: community ~ species")
print(t1)
cat("\n\ntest 2: effect of rhizosphere across all samples and seasons")
print(t2)
cat("\n\ntest 3: effect of rhizosphere per species")
print(t3)
cat("\n\ntest 4: change with season ")
print(t4)
cat("\n\ntest 5: change with season for every species")
print(t5)
sink()

#summary table for adonis results
table_adonis <-
  list(t1, t2, t3, t4, t5) %>%
    lapply(function(x) {
      x$aov.tab %>%
        as("data.frame") %>%
        cbind(formula  = x$call %>%
                as.character() %>% .[2]) %>%
        dplyr::select(formula, Df, `F.Model`, R2, `Pr(>F)`) %>%
        dplyr::mutate(strata = gsub(pattern = "^*\\$",
        replacement = "", x$call$strata)[3]) %>%
        .[1, ]
    }) %>%
  do.call(what = rbind) %>%
  dplyr::mutate(`F.Model` = round(`F.Model`, 2),
         R2 = round(R2, 2),
         `Pr(>F)` = round(`Pr(>F)`, 2))

flextable::flextable(table_adonis) %>%
  flextable::save_as_docx(path = "output/summary_table_adonis.docx")

write.csv(table_adonis, "output/adonis_summary_table.txt")
