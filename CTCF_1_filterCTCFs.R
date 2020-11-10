library(readr)
library(dplyr)
library(ggplot2)

### Read the database
ctcfs <- read_tsv("data/GSE85108/GSM2257816_CTCF_E0h_MCF7_peaks_hg38.bed", 
         col_names = c("chr", "start", "end", "id", "score"))

### How different are the peak sizes??
ctcfs <- ctcfs %>% mutate(size = end - start)

png(paste0("figures/ctcfs/size-density.png"), width = 600, height = 400)
ggplot(ctcfs, aes(x = size)) +
  scale_x_log10() +
  geom_density() + 
  theme_light()
dev.off()  

png(paste0("figures/ctcfs/size-boxplot.png"), width = 400, height = 600)
ggplot(ctcfs, aes(y = size)) +
  geom_boxplot() + 
  scale_y_log10() +
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

### Keep the 90% of the peaks, remove outliers
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

# Duplicate peaks
ctcfs %>% group_by(chr, start, end) %>% 
  tally() %>% filter(n > 1) %>% nrow()
# [1] 26
# Remove duplicate peaks, select highest score
ctcfs <- ctcfs %>% group_by(chr, start, end) %>% 
  filter(score == max(score))

# Remaining duplicates
ctcfs %>% group_by(chr, start, end) %>% 
  tally() %>% filter(n > 1) %>% nrow()
# [1] 6
ctcfs <- ctcfs %>% distinct(chr, start, end, .keep_all = T)

chrs <- as.character(c(seq(1:22), "X"))
ctcfs <- ctcfs %>% filter(chr %in% paste0("chr", chrs))

### Distance between peaks in chromosomes.
ctcfs_distance <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  ctcfs_chr <- ctcfs %>% filter(chr == paste0("chr", ch)) %>% 
    select(id, chr, start, end) %>%
    arrange(start) 
  ctcfs_chr$prev_end <- c(ctcfs_chr$start[1], ctcfs_chr$end[2:nrow(ctcfs_chr)-1])
  ctcfs_chr <- ctcfs_chr %>% mutate(distance = start - prev_end)
  return(ctcfs_chr)
})

ctcfs_distance <- plyr::ldply(ctcfs_distance)

### Around 100000 bps between peaks
png(paste0("figures/ctcfs/distance-density.png"), width = 400, height = 600)
ggplot(ctcfs_distance, aes(x = distance)) +
  geom_density() + 
  scale_x_log10() +
  theme_light()
dev.off()

png(paste0("figures/ctcfs/distance-boxplot-bychr.png"), width = 800, height = 400)
ggplot(ctcfs_distance, aes(y = distance, x = chr)) +
  geom_boxplot() + 
  scale_y_log10() +
  theme_light()
dev.off()

# Peak overlaps
peak_overlaps <- ctcfs_distance %>% filter(distance < 0) %>% 
  mutate(size = end - start, dif = size + distance)
nrow(peak_overlaps)
# [1] 23
# Remove them
ctcfs <- ctcfs %>% filter(!id %in% peak_overlaps$id) 

write_tsv(ctcfs, path = "data/ctcfs/filtered_ctcfs.tsv")
