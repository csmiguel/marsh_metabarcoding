###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel/marsh_metabarcoding
# May 2021
###.............................................................................
#GOAL: figure with rarefied sample seq depth
#PROJECT: marsh_metabarcoding
###.............................................................................
#
library(phyloseq)
library(vegan)
ps <- readRDS("data/intermediate/ps_t_noRep.rds")

pdf("output/rarefaction_curves.pdf", width = 5, height = 4, onefile = T)
vegan::rarecurve(otu_table(ps), step = 50, cex = 0.5,
          col = as.numeric(as.factor(sample_data(ps)$species)),
          lwd = 2, ylab = "ASVs")
legend("bottomright", legend = levels(as.factor(sample_data(ps)$species)),
  col = 1:nlevels(as.factor(sample_data(ps)$species)),
  lty = 1, lwd = 2, bty = "n", cex = 0.6)
dev.off()
