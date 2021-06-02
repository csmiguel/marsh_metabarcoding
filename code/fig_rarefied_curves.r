# figure of rarefied samples
library(vegan)
ps_marsh <- readRDS("data/intermediate/ps_marsh_filt.rds")
ps_greenh <- readRDS("data/intermediate/ps_greenh_filt.rds")

layout_matrix_1 <- matrix(1:2, ncol = 1)
pdf("output/rarefaction_curves.pdf", width = 5, height = 8, onefile = T)
layout(layout_matrix_1)
rarecurve(otu_table(ps_marsh), step = 50, cex = 0.5,
          col = as.numeric(as.factor(sample_data(ps_marsh)$species)),
          lwd = 2, ylab = "ASVs")
rarecurve(otu_table(ps_greenh), step = 50, cex = 0.5,
          col = as.numeric(as.factor(sample_data(ps_greenh)$species)),
          lwd = 2, ylab = "ASVs")
dev.off()
