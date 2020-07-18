library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)

conds <- c("healthy", "luma")
label_conds <- c("Healthy", "LumA")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
names(colors) <- conds

all_assort_vals <- lapply(conds, function(cond){
  assort_vals <- read_tsv(paste0("data/assortativity/", cond, "-chr-assortativity.tsv"))
  assort_vals$cond <- cond
  return(assort_vals)
})
all_assort_vals <- bind_rows(all_assort_vals)

# Plot
p <- all_assort_vals %>%
  ggplot( aes(x=cond, y=totalfrac, fill=cond, color=cond)) +
  geom_violin(width=1.5, size=0.2) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_base() +
  theme(
    legend.position="none"
  ) +
  scale_x_discrete(labels = label_conds)+
  xlab("") +
  ylab("Fraction of intra links")

png(paste0("figures/communities/community-fraction-intra-links.png"), width = 1200, height = 800)
p
dev.off() 

luma_exp <- read_tsv(paste0("data/assortativity/luma-exp-assortativity.tsv"))

p <-ggplot(luma_exp, aes(x = totalfrac, fill="luma")) + 
  geom_histogram(bins = 50) +
  theme_base() +
  theme(
    legend.position="none"
  ) +
  scale_fill_manual(values = colors["luma"]) +
  xlab("Fraction of same expression links") +
  ylab("Frequency")

png(paste0("figures/communities/community-fraction-exp-links.png"), width = 1200, height = 800)
p
dev.off() 