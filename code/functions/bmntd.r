
# function to shuffle otu_table
tax_shuffle <- function(x) {
  # x, data.frame with cols being taxa
  colnames(x) <- sample(colnames(x))
  return(x)
}

# function to carry out pairwise computing of bMNTD permutations
permutate_bmntd <- function(physeq = NULL, nperm = nperm) {
  # create input for picante
  otu_t <- as.data.frame(otu_table(physeq))
  ctree <- cophenetic(phy_tree(physeq))
  # check
  stopifnot(nsamples(physeq) == 2)
  # run permutations
  temp <-
    lapply(1:nperm, function(x) {
      sotu_t <- tax_shuffle(otu_t)
      tempdist <-
        picante::comdistnt(comm = sotu_t,
                           dis = ctree,
                           abundance.weighted = T,
                           exclude.conspecifics = T)
      melt(as.matrix(tempdist), varnames = c("row", "col")) %>% .[3,]
    }) %>%
    do.call(what = "rbind")
  return(temp)
}
