###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# April 2021
###.............................................................................
#GOAL: to produce a table with reads retained at every step of the
# filterning process and ASV determination
#PROJECT: spartina-metarizo
###.............................................................................
library(dada2)
library(dplyr)

#filtered and trimmed, in and out
truncated <- readRDS("data/intermediate/truncated-reads.rds") %>%
            `rownames<-`(sapply(strsplit(rownames(.), "_"), `[`, 1))
#dada asvs
load("data/intermediate/dada.Rdata")
#merged seqs
load("data/intermediate/mergers.Rdata")
#merged seqs size filtered
load("data/intermediate/seqtab2.Rdata")
#OTU tables: no chimeras
load("data/intermediate/seqtabNoC.Rdata")

# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
names_track <- c("sample", "input", "filtered", "denoised",
  "merged", "tabled", "nonchim", "variants")

#marsh
track_marsh <-
        truncated[grepl(pattern = "S-[0-9]|BPCR", rownames(truncated)), ] %>%
        as.data.frame() %>% rename(reads_in = 1, reads_out = 2) %>%
        tibble::rownames_to_column("sample") %>%
        dplyr::left_join(
          sapply(marsh_dadaFs, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          sapply(mergers_marsh, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtab_marsh_2) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtabNoC_marsh) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          apply(seqtabNoC_marsh, 1, function(x) {sum(x > 0)}) %>%
          as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        setNames(names_track)

#Spartina maritima
track_greenh <-
        truncated[grepl(pattern = "MA-[0-9]|BPCR", rownames(truncated)), ] %>%
        as.data.frame() %>% rename(reads_in = 1, reads_out = 2) %>%
        tibble::rownames_to_column("sample") %>%
        dplyr::left_join(
          sapply(greenh_dadaFs, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          sapply(mergers_greenh, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtab_greenh_2) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtabNoC_greenh) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          apply(seqtabNoC_greenh, 1, function(x) {sum(x > 0)}) %>%
          as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        setNames(names_track)

#write results
write.csv(track_greenh, "output/tracked_reads_greenhouse_exp.csv")
write.csv(track_marsh, "output/tracked_reads_marsh-plants_exp.csv")
#write to excel
xlsx::write.xlsx(track_greenh,
  file = "output/results.xls",
  sheetName = "tracked_reads_greenhouse_exp")
xlsx::write.xlsx(track_marsh,
  file = "output/results.xls",
  sheetName = "tracked_reads_marsh", append = T)
