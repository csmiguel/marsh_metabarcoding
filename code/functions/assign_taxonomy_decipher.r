#assign taxonomy using DECIPHER
# a more detailed description in http://benjjneb.github.io/dada2/tutorial.html
# DECIPHER::idTaxa uses a better performing
# classification algorithm https://doi.org/10.1186/s40168-018-0521-5.

assign_taxonomy_decipher <- function(seqtabdada2 = NULL,
  training_set = NULL,
  nprocessors = 1,
  in_parallel = F) {
    #seqtabdada2, is a OTU like table produced by makeSequenceTable
    # (and possibly filtered, eg removeBimeraDenovo)
    #training_set, is an R object formatte for DECIPHER. Downloaded from:
    # "http://www2.decipher.codes/Classification/TrainingSets"
    #nprocessors, number of processors to use. In MAC the max is 1
    # see more info here: https://github.com/benjjneb/dada2/issues/1333
    #required libraries
    require(dada2)
    require(DECIPHER)
    require(dplyr)
    #assertions
    assertthat::assert_that(
      all(class(training_set) == c("Taxa", "Train")),
      msg = "check format of training_set"
    )
    assertthat::assert_that(
      all(class(seqtabdada2) == c("matrix", "array")),
      msg = "check format of seqtabdada2"
    )
    # extract sequences from OTU object
    dna <-
      Biostrings::DNAStringSet(dada2::getSequences(seqtabdada2))
    if (in_parallel == F) {
    # DNAStringSet from the ASVs
    ids <- DECIPHER::IdTaxa(dna,
                            training_set,
                            #for MACs, only 1 processor (see below)
                            processors = nprocessors,
                            strand = "top",
                            verbose = T)
    } else if (in_parallel == T) {
      #do parallel processing
      assertthat::assert_that(nprocessors > 1,
      msg = "If 1 processor is selected, use in_parallel = F")
      require(parallel)
      cat("Running in parallele\n", nprocessors, "created")
      #function to split dna object into ~similar size chunks
      # for parallel processing
      split_vector <- function(for_n_clusters = NULL,
        query_seqs = NULL) {
        h <-
          rep(letters[1:for_n_clusters],
            length.out = length(query_seqs),
            each = ceiling(length(query_seqs) /
                    for_n_clusters))[1:length(query_seqs)]
        as.factor(h)
      }
      #create clusters
      cl <- parallel::makeCluster(nprocessors,
         #https://github.com/rstudio/rstudio/issues/6692
                        setup_strategy = "sequential")
      #add desired environment variables to clusters
      parallel::clusterExport(cl, c("dna", "training_set"), envir = environment())
      #list with positions of dna seqs from the object dna to ID in each cluster
      hh <-
        seq_along(dna) %>%
        split(f = split_vector(for_n_clusters = nprocessors, query_seqs = dna))
      # run IdTaxa in parallel
      ids <-
        parallel::parLapply(cl, hh, function(x) {
          DECIPHER::IdTaxa(dna[x],
                                  training_set,
                                  processors = 1,
                                  strand = "top",
                                  verbose = T)
        }) %>% do.call(what = c) %>% setNames(NULL)
      parallel::stopCluster(cl)
          }
    # ranks of interest
    ranks <- c("domain", "phylum", "class", "order",
               "family", "genus", "species")
    # Convert the output object of class "Taxa" to a matrix analogous to
    # the output from assignTaxonomy
    taxid <-
      ids %>%
      sapply(function(x) {
      m <- match(ranks, x$rank)
      taxa <- x$taxon[m]
      taxa[startsWith(taxa, "unclassified_")] <- NA
      taxa
    }) %>% t()
    colnames(taxid) <- ranks
    rownames(taxid) <- dada2::getSequences(seqtabdada2)
    taxid
  }
