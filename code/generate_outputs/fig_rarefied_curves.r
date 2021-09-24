# figure of rarefied samples
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
