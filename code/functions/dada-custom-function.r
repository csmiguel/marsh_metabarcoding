#define function for facilitate the processing of data from both experiments:
# 1. marsh: all 5 species in the marsh soil.
# 2. greenh: greenhouse experiment with Spartina maritima.
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
