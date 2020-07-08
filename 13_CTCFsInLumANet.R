library(readr)
library(dplyr)
library(ggplot2)
library(IRanges)

ctcfs <- read_tsv("data/class_ctcfs.tsv", col_types = cols_only(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  id = col_character(),
  score = col_double()))

genes <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(
                            ensemblID = col_character(),
                            chr = col_character(),
                            start = col_integer(),
                            end = col_integer()
                          ))
genes <- genes %>% dplyr::rename(id = ensemblID)

promoters <- genes %>% mutate(end = start + 5000, start = start - 5000)
extended_regions <- genes %>% mutate(start = start - 5000, end = end)
gene_bodies <- genes %>% mutate(start = start + 5000) %>% 
  filter(end > start)

chrs <- as.character(c(seq(1:22), "X"))

luma_ctfs <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  chr_ctcfs <- ctcfs %>% filter(chr == ch)
  chr_gene_bodies <- gene_bodies %>% filter(chr == ch)
  chr_promoters <- promoters %>% filter(chr == ch)
  chr_extended_regions <- extended_regions %>% filter(chr == ch)
  
  cranges <- IRanges(start = chr_ctcfs$start, end = chr_ctcfs$end, names = chr_ctcfs$id)
  granges <- IRanges(start = chr_gene_bodies$start, end = chr_gene_bodies$end, names = chr_gene_bodies$id)
  pranges <- IRanges(start = chr_promoters$start, end = chr_promoters$end, names = chr_promoters$id)
  eranges <- IRanges(start = chr_extended_regions$start, end = chr_extended_regions$end, names = chr_extended_regions$id)
  
  #ctcfs in promotoers and gene bodies
  chr_ctcfs <- chr_ctcfs %>% mutate(phits = countOverlaps(cranges, pranges, type = "within"), 
                                    ghits = countOverlaps(cranges, granges, type = "within"))
                         
  phits <- as.data.frame(findOverlaps(cranges, pranges, type = "within")) %>% 
    mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_promoters$id[subjectHits], type = "promoter")
  
  ghits <- as.data.frame(findOverlaps(cranges, granges, type = "within")) %>% 
    mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_gene_bodies$id[subjectHits], type = "gene")
  
  hits <- bind_rows(phits, ghits) %>% select(ctcf, sequence, type)
  
  # Distance to the extended region
  dist <- distanceToNearest(cranges, eranges)
  chr_ctcfs$near_region <- chr_extended_regions$id[to(dist)]
  chr_ctcfs$near_distance <- unlist(elementMetadata(dist))
  return(list(ctcfs = chr_ctcfs, hits = hits))
})

class_ctcfs <- plyr::ldply(lapply(luma_ctfs, "[[", "ctcfs"))
ctcfs_hits <- plyr::ldply(lapply(luma_ctfs, "[[", "hits"))

class_ctcfs <- class_ctcfs %>% mutate(type = case_when(phits > 0 ~ "promoter",
                                                   ghits > 0 ~ "gene",
                                                   TRUE ~ "intergen"))

class_ctcfs %>% select(type) %>% table()
# gene intergen promoter 
# 1749    18184      420 


png(paste0("figures/ctcfs/class-barplot-luma.png"), width = 600, height = 800)
ggplot(class_ctcfs, aes(x = type, fill = type)) +
  geom_bar() +
  theme_bw() +
  scale_fill_viridis_d() 
dev.off()

ctcfs_count <- bind_rows(
  class_ctcfs %>% filter(near_distance <= 100000) %>% select(type) %>% mutate(distance = "100k"),
  class_ctcfs %>% filter(near_distance <= 75000) %>% select(type) %>% mutate(distance = "75k"),
  class_ctcfs %>% filter(near_distance <= 50000) %>% select(type) %>% mutate(distance = "50k"),
  class_ctcfs %>% filter(near_distance <= 25000) %>% select(type) %>% mutate(distance = "25k"),
  class_ctcfs %>% filter(near_distance <= 15000) %>% select(type) %>% mutate(distance = "15k"),
  class_ctcfs %>% filter(near_distance <= 10000) %>% select(type) %>% mutate(distance = "10k"),
  class_ctcfs %>% filter(near_distance <= 5000) %>% select(type) %>% mutate(distance = "5k"),
  class_ctcfs %>% filter(near_distance <= 1000) %>% select(type) %>% mutate(distance = "1k"))
ctcfs_count <- ctcfs_count %>% group_by(type, distance) %>% tally()
ctcfs_count$distance <- factor(ctcfs_count$distance, 
                                  levels = c("1k", "5k", "10k", "15k", "25k", "50k", "75k", "100k"))

png(paste0("figures/ctcfs/class-barplot-luma-distance.png"), width = 800, height = 800)
ggplot(ctcfs_count, aes(x = distance,  y = n, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() 
dev.off()

## digamos que 25k 
class_ctcfs_25 <- class_ctcfs %>% filter(near_distance <= 25000)

png(paste0("figures/ctcfs/class-barplot-luma-distance-bychr.png"), width = 1200, height = 1200)
ggplot(class_ctcfs_25, aes(x = type, fill = type)) +
  geom_bar() +
  theme_bw() +
  scale_fill_viridis_d() + 
  facet_wrap(~chr)
dev.off()

write_tsv(class_ctcfs, path = "data/ctcfs_in_luma.tsv")
