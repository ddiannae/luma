library(readr)
library(dplyr)
library(ggplot2)
library(IRanges)

### Read the classified CTCFs
ctcfs <- read_tsv("data/class_ctcfs_5000_500.tsv", col_types = cols(
  chr = col_character()))

ctcfs %>% select(type) %>% table()
# -5000, 500
# gen  intergen    other promoter 
# 8047     8979     2459      868
# intergen = 11438
# 
# -1000, 100
# gen intergen    other promoter 
# 8270     9393     2523      167 
# intergen = 11916
ctcfs_hits <- read_tsv("data/ctcfs_hits_5000_500.tsv")
genes <- read_tsv("data/luma-intra-vertices.tsv", 
                          col_types = cols_only(
                            ensemblID = col_character(),
                            chr = col_character(),
                            start = col_integer(),
                            end = col_integer()
                          ))
genes <- genes %>% dplyr::rename(id = ensemblID)

luma_ctcfs <- ctcfs_hits %>% 
  semi_join(genes, by = c("sequence" = "id")) %>% 
  select(ctcf) %>% unlist(use.names = FALSE) %>% unique()

## Remove ctcfs that lie on gene promoters or gene bodies in genes that are not 
## in our network. Keep the ones in the intergenic region of genes to check if they include
## other peaks at a certain distance
ctcfs <- bind_rows(
  ctcfs %>% 
    filter(id %in% luma_ctcfs),
  ctcfs %>% 
  filter(!id %in% luma_ctcfs, type  %in% c("intergen", "other"))
  )

promoters <- genes %>% mutate(end = start + 500, start = start - 5000)
extended_regions <- genes %>% mutate(start = start - 5000, end = end)
gene_bodies <- genes %>% mutate(start = start + 500) %>% filter(end > start)

chrs <- as.character(c(seq(1:22), "X"))

### Reclasify
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
                         
  phits <- as.data.frame(findOverlaps(cranges, pranges, type = "within")) 
  if(nrow(phits) > 0) {
    phits <- phits  %>% 
      mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_promoters$id[subjectHits], type = "promoter")
  }
 
  ghits <- as.data.frame(findOverlaps(cranges, granges, type = "within"))
  if(nrow(ghits) > 0) {
    ghits <- ghits %>% 
      mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_gene_bodies$id[subjectHits], type = "gene")
  }
   
  hits <- bind_rows(phits, ghits)
  
  if(nrow(hits) > 0) {
    hits <- hits %>% select(ctcf, sequence, type)
  }
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
# -5000, 500
# gene intergen promoter 
# 1343    11438      177
# 
# -1000, 100
# gene intergen promoter 
# 1370    11916       48 

png(paste0("figures/ctcfs/class-barplot-luma-5000-500.png"), width = 600, height = 800)
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

png(paste0("figures/ctcfs/class-barplot-luma-distance-5000-500.png"), width = 800, height = 800)
ggplot(ctcfs_count, aes(x = distance,  y = n, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() 
dev.off()

## digamos que 50k 
class_ctcfs_50 <- class_ctcfs %>% filter(near_distance <= 50000)

class_ctcfs_50 %>% select(type) %>% table()
# -5000, 500
# gene intergen promoter 
# 1343      887      177 

# -1000, 100
# gene intergen promoter 
# 1370     1027       48 
plot_ctcfs <- class_ctcfs_50 %>% select(type, chr) %>% 
  bind_rows(genes %>% select(chr) %>% mutate(type = "mrna")) %>% 
  mutate(type = ordered(as.factor(type), levels = c("mrna", "gene", "promoter", "intergen")),
         chr = ordered(chr, levels = chrs))
  
png(paste0("figures/ctcfs/ctcfs-by-chr-type-5000-500.png"), width = 1200, height = 1200)
ggplot(plot_ctcfs, aes(x = type, fill = type)) +
  geom_bar() +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank()) +
  scale_fill_viridis_d(option = "C", name = "Type",
                       labels = c("Genes", "CTCF bs in genes", "CTCF bs in promoters", "CTCF bs in intergenic\n region" )) + 
  facet_wrap(~chr) +
  xlab("") +
  ylab("") +
  ggtitle("Genes and CTCF bs in intra-chromosomal communities") 
dev.off()

write_tsv(class_ctcfs, path = "data/ctcfs_in_intra_luma_5000_500.tsv")
write_tsv(ctcfs_hits, path  = "data/ctcfs_hits_in_intra_luma_5000_500.tsv")
write_tsv(class_ctcfs_50, path  = "data/ctcfs_in_intra_luma_50k_5000_500.tsv")

