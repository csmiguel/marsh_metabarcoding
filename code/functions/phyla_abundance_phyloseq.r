plot_abundances <- function(phylo) {
  require(ggplot2)
  prev <-
    apply(X = otu_table(phylo),
         MARGIN = ifelse(taxa_are_rows(phylo), yes = 1, no = 2),
         FUN = function(x){sum(x > 0)})
  # add taxonomy and total read counts to this data.frame
  prev <- data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(phylo),
                    tax_table(phylo))

  # get prevalence per phylum
  prevalence_per_phylum <-
    plyr::ddply(prev, "phylum", function(df1){
      cbind(mean(df1$Prevalence), sum(df1$Prevalence))
      })
  # plot abundances vs prevalence
  plot1 <-
    prev %>%
    ggplot(aes(TotalAbundance, Prevalence / nsamples(phylo), color = phylum)) +
    # Include a guess for parameter
    geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +
    xlab("Total Abundance") +
    ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~phylum) +
    theme_minimal() +
    theme(legend.position = "none")
  return(list(prev, prevalence_per_phylum, plot1))
  }
