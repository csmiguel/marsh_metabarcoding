# figure of rarefied samples
library(vegan)
ps <- readRDS("data/intermediate/ps_filt.rds")

pdf("output/rarefaction_curves.pdf", width = 5, height = 8, onefile = T)
rarecurve(otu_table(ps), step = 50, cex = 0.5,
          col = as.numeric(as.factor(sample_data(ps)$species)),
          lwd = 2, ylab = "ASVs")
dev.off()
