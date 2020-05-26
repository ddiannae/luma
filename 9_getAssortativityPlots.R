library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)

conds <- c("healthy", "luma")

all_assort_vals <- lapply(conds, function(cond){
  assort_vals <- read_tsv(paste0("data/assortativity/", cond, "-chr-assortativity.tsv"))
  assort_vals$cond <- cond
  return(assort_vals)
})
all_assort_vals <- bind_rows(all_assort_vals)

# Plot
p <- all_assort_vals %>%
  ggplot( aes(x=cond, y=totalfrac, fill=cond, color=cond)) +
  geom_violin(width=2.1, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_base() +
  theme(
    legend.position="none"
  ) +
 xlab("") +
 ylab("Fraction of intra links")

luma_exp <- read_tsv(paste0("data/assortativity/luma-exp-assortativity.tsv"))

p <-ggplot(luma_exp, aes(x = totalfrac)) + 
  geom_histogram() +
  theme_base() +
  scale_fill_viridis(discrete=TRUE) +
  xlab("Fraction of same category links") +
  ggtitle("Distribution of expression assortativity") +
  ylab("Frequency")
