###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# May 2021
###.............................................................................
#GOAL: compute alpha diversity bootstrapped
#PROJECT: spartina-metarizo
###.............................................................................
library(phyloseq)
library(ggplot2)
library(dplyr)

#read filtered phyloseq
ps <- readRDS("data/intermediate/ps_marsh_filt.rds")
