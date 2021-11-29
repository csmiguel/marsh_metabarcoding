###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez AT gmail DOT com // miguelcamachosanchez.weebly.com
# https://scholar.google.com/citations?user=1M02-S4AAAAJ
# https://orcid.org/0000-0002-6385-7963
# github.com/csmiguel
# August 2021
###.............................................................................
#GOAL: guess functions from ASVs
#PROJECT: spartina-metarizo
###...........................................................................

library(Tax4Fun2)
#  Tax4fun2 is compatible with the Silva database. Otherwise, Picrust2 is
# compatible with the Greengenes database.
Tax4Fun2::buildReferenceData(install_suggested_packages	= F)
#https://github.com/bwemheu/Tax4Fun2
#https://www.genome.jp/kegg/pathway.html
