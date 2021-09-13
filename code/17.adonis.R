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

# variable with permutations
nperm <- 9999

# read phyloseq
ps <- readRDS("data/intermediate/ps_log.rds")

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
