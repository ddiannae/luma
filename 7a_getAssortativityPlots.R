library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)

conds <- c("healthy", "luma")
label_conds <- c("Healthy", "Luminal A")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
names(colors) <- conds

labels_alg <- c("Fast Greedy", "Infomap", "Leading Eigenvector", "Louvain")
algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector","multi_level")

names(labels_alg) <- algorithms

m <- lapply(algorithms, function(algrthm) {
  
  all_assort_vals <- lapply(conds, function(cond){
    assort_vals <- read_tsv(paste0("data/assortativity/", cond, 
                                   "-chr-assortativity-", algrthm,".tsv"))
    assort_vals$cond <- cond
    return(assort_vals)
  })
  
  all_assort_vals <- bind_rows(all_assort_vals)
  
  # Plot
  p <- all_assort_vals %>%
    ggplot( aes(x=cond, y=totalfrac, fill=cond, color=cond)) +
    geom_violin() +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_base() +
    theme(
      legend.position="none"
    ) +
    scale_x_discrete(labels = label_conds)+
    xlab("") +
    ylab("Fraction of intra-chromosomal links") +
    theme(text = element_text(size = 20), axis.text.x = element_text(size = 28),
          axis.title.y = element_text(size = 25), plot.background=element_blank())
  
  png(paste0("figures/assortativity/community-fraction-intra-links-", algrthm, ".png"), 
      width = 1000, height = 500)
  print(p)
  dev.off()

  # Plot
  p <- all_assort_vals %>%
    ggplot( aes(x=cond, y=diffraction, fill=cond, color=cond)) +
    geom_violin() +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_base() +
    theme(
      legend.position="none"
    ) +
    scale_x_discrete(labels = label_conds)+
    xlab("") +
    ylab("Chromosomal Assortativity") +
    theme(text = element_text(size = 20), axis.text.x = element_text(size = 28),
          axis.title.y = element_text(size = 25), plot.background=element_blank())

  png(paste0("figures/assortativity/community-diff-intra-links-", algrthm, ".png"),
      width = 1000, height = 500)
  print(p)
  dev.off()
  
  all_assort_vals$algrthm <-labels_alg[algrthm]
  return(all_assort_vals)
})

m <- bind_rows(m)

p <- m %>%
  ggplot(aes(x=cond, y=diffraction, fill=cond, color=cond)) +
  geom_violin() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_base() +
  theme(
    legend.position="none"
  ) +
  scale_x_discrete(labels = label_conds)+
  xlab("") +
  ylab("Chromosomal Assortativity") +
  theme(text = element_text(size = 20), axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 16), plot.background=element_blank()) +
  facet_wrap(~algrthm, ncol = 1, labeller = label_value) 

png(paste0("figures/assortativity/community-diff-intra-links.png"),
    units="in", width=5, height=8, res=300)
print(p)
dev.off()

luma_exp <- read_tsv(paste0("data/assortativity/luma-exp-assortativity.tsv"))

p <-ggplot(luma_exp, aes(x = totalfrac, fill="luma")) + 
  geom_histogram(bins = 50) +
  theme_base() +
  theme(
    legend.position="none", text = element_text(size = 16)
  ) +
  scale_fill_manual(values = colors["luma"]) +
  xlab("Fraction of same expression links") +
  ylab("Frequency")

png(paste0("figures/assortativity/community-fraction-exp-links.png"), 
    units="in", width=5, height=8, res=300)
print(p)
dev.off() 
