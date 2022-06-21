library(ggplot2)
library(phyloseq)
library(dplyr)
library(reshape2)

# import results
ft <-
  read.delim("output/functional_table.tsv")
# filter rows with data
# T/F vector
filt_vector <- ft %>% select(-group) %>% rowSums() > 0
# filter rows and melt
df1 <-
  ft[filt_vector, ] %>%
  reshape2::melt(value.name = "freq") %>%
  mutate(sample_name = gsub("\\.", "-", variable)) %>%
  select(-variable)
# metadata from samples
sdata <-
  sample_data(ps) %>%
  as("data.frame") %>%
  select(-rep) %>%
  arrange(species, season, rhizosphere)
# desired order of samples
order_samples <- sdata$sample_name
match(plotdata$sample_name, order_samples)
factor(mydf$task, levels = c("up", "down", "left", "right", "front", "back"))
plotdata <-
  left_join(df1, sdata, by = "sample_name") %>%
  mutate(sample_name = factor(sample_name,
                              levels = order_samples))

bubble_plot <-
  ggplot(plotdata, aes(y = group, x = sample_name)) +
  geom_point(aes(size=ifelse(freq == 0, NA, freq), fill = species),
             shape=21, color="black",
             alpha = as.numeric(as.factor(plotdata$rhizosphere))/2) +
  guides(size = guide_legend(title = "relative frequency")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

ggsave("output/faprotax.pdf", bubble_plot, height = 7, width = 9)
