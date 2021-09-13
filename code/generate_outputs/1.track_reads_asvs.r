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
truncated <-
  readRDS("data/intermediate/truncated_reads.rds") %>%
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
names_track <-
  c("sample", "input", "filtered", "denoised",
  "merged", "tabled", "nonchim", "variants")

track <-
        truncated[grepl(pattern = "S-[0-9]|BPCR", rownames(truncated)), ] %>%
        as.data.frame() %>% rename(reads_in = 1, reads_out = 2) %>%
        tibble::rownames_to_column("sample") %>%
        dplyr::left_join(
          sapply(dadaFs, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          sapply(mergers, getN) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtab2) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          rowSums(seqtabNoC) %>% as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        dplyr::left_join(
          apply(seqtabNoC, 1, function(x) {sum(x > 0)}) %>%
          as.data.frame() %>% rename(x = 1) %>%
            tibble::rownames_to_column("sample"), by = "sample") %>%
        setNames(names_track)


#write results
write.csv(track, "output/tracked_reads.csv")

#write to excel
xlsx::write.xlsx(track,
  file = "output/results.xls",
  sheetName = "tracked_reads", append = T)
