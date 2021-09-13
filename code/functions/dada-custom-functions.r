#define function for facilitate the processing of data from both experiments
# so they dropped in different pools for ASV determination.
########
dada_1 <- function(dereplicated = NULL, pattern_derep = NULL, errors = NULL) {
  # dereplicated is a dada de-replicated object.
  # pattern is the pattern to grep from dereplicated names.
  # error is a learnt errors object.
  grep(x = names(dereplicated), pattern = pattern_derep) %>%
  dereplicated[.] %>%
    dada2::dada(
      err = errors,
      multithread = 2,
      pool = T)
    }
##########
mergerPairs_1 <- function(dada_f = NULL,
                          dada_r = NULL,
                          derep_f = derepFs,
                          derep_r = derepRs,
                          pattern_sampleset = NULL) {
  # mergePairs
  # pattern is the pattern to grep from dereplicated names.
  assertthat::assert_that(class(dada_f) == "list" &
                          class(dada_r) == "list" &
                          class(derep_f) == "list" &
                          class(derep_r) == "list")
# filter derep samples that match pattern
  derep_f_filt <-
        grep(x = names(derep_f),
             pattern = pattern_sampleset) %>%
        derep_f[.]
  derep_r_filt <-
        grep(x = names(derep_r),
             pattern = pattern_sampleset) %>%
        derep_r[.]
# run mergePairs
  dada2::mergePairs(
        dada_f, derep_f_filt,
        dada_r, derep_r_filt,
        verbose = TRUE,
        minOverlap = 20,
        maxMismatch = 0)
    }
