###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2021
###.............................................................................
#GOAL: assign taxonomy to ASVs using Silva v138 database:
#PROJECT: spartina-metarizo
###.............................................................................
library(dada2)
library(dplyr)
library(DECIPHER)

#load OTU tables
load("data/intermediate/seqtabNoC.Rdata")
#custom function to assign taxonomy based on IdTaxa
source("code/functions/assign_taxonomy_decipher.r")

#a formatted version of the latest release of SILVA is maintained for dada2 in
#Zenodo 10.5281/zenodo.4587955. However, from dada2, Benjamin recomends using
#DECIPHER curated database and the function DECIPHER::idTaxa instead of
#dada2::assignTaxonomy. After, I tested both, idTaxa does not use so much RAM.

#1. get training set for SILVA.
#paths to training set from SILVA
db_silva <- "SILVA_SSU_r138_2019.RData"
dest_path <- "data/raw"

#download training dataset SILVA formatted for DECIPHER if not yet downloaded
url_silva <- "http://www2.decipher.codes/Classification/TrainingSets"
if (!file.exists(file.path(dest_path, db_silva))) {
  options(timeout = 1000) #increase the default max time allowed to download a
  download.file(file.path(url_silva, db_silva),
                destfile = file.path(dest_path, db_silva))
}
#load SILVA trainingSet
load(file.path(dest_path, db_silva))

#2. assign taxonomy
# marsh
taxid_marsh <-
  assign_taxonomy_decipher(seqtabNoC_marsh,
  training_set = trainingSet,
  nprocessors = 1,
  in_parallel = F)
# greenhouse
taxid_greenh <-
  assign_taxonomy_decipher(seqtabNoC_greenh,
  training_set = trainingSet,
  nprocessors = 1,
  in_parallel = F)
#taxid object has the below structure. rownames are seqs and taxonomy in cols.
# (rownames) domain  phylum    class      order     family    genus     species
# ACCTAT… Bacter… Campilobac… Campylob… Campylo… Sulfurovace… Sulfurovum    NA
# ACCTCT… Bacter… Proteobact… Alphapro… Defluvi… Defluviicoc… Defluviicocc… NA
# ACCTGT… Bacter… Desulfobac… Desulfob… Desulfo… Desulfosarc… NA            NA

#save seqs with assigned taxonomy
saveRDS(taxid_marsh, file = "data/intermediate/ids_marsh.rds")
saveRDS(taxid_greenh, file = "data/intermediate/ids_greenh.rds")
