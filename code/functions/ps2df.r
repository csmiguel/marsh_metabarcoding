ps2df <- function(ps = NULL) {
  #convert phyloseq slots to list of dataframes
  #converts OTU Table, Sample Data and Taxonomy Table, dataframes
  #output is a list
  require(phyloseq)
  require(dplyr)
  assertthat::assert_that(
    all(c("otu_table", "tax_table", "sam_data") %in% getslots.phyloseq(ps)),
    msg = "Missing slots in phyloseq")
  #otu table
 ps_otu_table <-
   otu_table(ps) %>%
   as.data.frame()
  #sample data
  h <- sample_data(ps)
  ps_sample_data <-
    h@.Data %>%
    sapply(function(x) x) %>%
    as.data.frame() %>%
    setNames(names(h)) %>%
    `rownames<-`(rownames(h))
  #taxonomy table
  ps_tax_table <-
    tax_table(ps) %>%
    as.data.frame()
  list(otu_table = ps_otu_table,
       sample_data = ps_sample_data,
       tax_table = ps_tax_table)
  }
