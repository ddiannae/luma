library(readr)
library(dplyr)
library(ggplot2)

ctcfs <- read_tsv("data/GSE85108/GSM2257816_CTCF_E0h_MCF7_peaks_hg38.bed", 
         col_names = c("chr", "start", "end", "id", "score"))

ctcfs <- ctcfs %>% mutate(size = end - start)

png(paste0("figures/ctcfs/size-density.png"), width = 600, height = 400)
ggplot(ctcfs, aes(x = size)) +
  geom_density() + 
  theme_light()
dev.off()  

png(paste0("figures/ctcfs/size-boxplot.png"), width = 400, height = 600)
ggplot(ctcfs, aes(y = size)) +
  geom_boxplot() + 
  theme_light()
dev.off()  

mean(ctcfs$size)
#[1] 346.8707
quantile(ctcfs$size)
#0%   25%   50%   75%  100% 
#131   245   312   379 54545
quantile(x = ctcfs$size, probs = c(0.9, 0.95, 0.98))
# 90%    95%    98% 
# 482.00 577.35 773.14 
ctcfs <- ctcfs %>% filter(size <= 482)

png(paste0("figures/ctcfs/size-density-filtered.png"), width = 600, height = 400)
ggplot(ctcfs, aes(x = size)) +
  geom_density() + 
  theme_light()
dev.off()  

png(paste0("figures/ctcfs/size-boxplot-filtered.png"), width = 400, height = 600)
ggplot(ctcfs, aes(y = size)) +
  geom_boxplot() + 
  theme_light()
dev.off()  

png(paste0("figures/ctcfs/score-density.png"), width = 600, height = 400)
ggplot(ctcfs, aes(x = score)) +
  geom_density() + 
  theme_light()
dev.off()  

png(paste0("figures/ctcfs/score-boxplot.png"), width = 400, height = 600)
ggplot(ctcfs, aes(y = score)) +
  geom_boxplot() + 
  theme_light()
dev.off() 

write_tsv(ctcfs, path = "data/filtered_ctcfs.tsv")
