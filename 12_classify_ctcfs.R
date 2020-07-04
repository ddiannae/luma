library(readr)
library(dplyr)
library(plyr)
library(stringr)
library(IRanges)

chrs <- as.character(c(seq(1:22), "X"))
annot <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt",
                  col_names = c("id", "chr", "start", "end", "gc", "type", "symbol"),
                  col_types = cols(col_character(), col_character(), col_integer(), 
                                   col_integer(), col_double(), col_character(), col_character()),
                  skip = 1)
annot <- annot %>% filter(chr %in% chrs)
protein_coding <- annot %>% filter(type == "protein_coding")
non_protein_coding <- annot %>% filter(type != "protein_coding")

ctcfs <- read_tsv("data/filtered_ctcfs.tsv")
ctcfs <- ctcfs %>% mutate(chr = str_remove(chr, "chr")) %>% 
  filter(chr %in% chrs)

promoters <- protein_coding %>% mutate(end = start + 5000, start = start - 5000)
genes <- protein_coding %>% mutate(start = start + 5000)
genes <- genes %>% filter(end > start)

ctcfs_hits <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  chr_ctcfs <- ctcfs %>% filter(chr == ch)
  chr_genes <- genes %>% filter(chr == ch)
  chr_promoters <- promoters %>% filter(chr == ch)
  chr_other <- non_protein_coding %>% filter(chr == ch)
  
  cranges <- IRanges(start = chr_ctcfs$start, end = chr_ctcfs$end, names = chr_ctcfs$id)
  granges <- IRanges(start = chr_genes$start, end = chr_genes$end, names = chr_genes$id)
  pranges <- IRanges(start = chr_promoters$start, end = chr_promoters$end, names = chr_promoters$id)
  oranges <- IRanges(start = chr_other$start, end = chr_other$end, names = chr_other$id)
  
  chr_ctcfs <- chr_ctcfs %>% mutate(phits = countOverlaps(cranges, pranges, type = "within"), 
                            ghits = countOverlaps(cranges, granges, type = "within"),
                            ohits = countOverlaps(cranges, oranges, type = "within"))
  
  promoters <- as.data.frame(findOverlaps(cranges, pranges, type = "within"))
  promoters <- promoters %>% 
    mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_promoters$id[subjectHits], type = "promoter")
  
  genes <- as.data.frame(findOverlaps(cranges, granges, type = "within"))
  genes <- genes %>% 
    mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_genes$id[subjectHits], type = "gene")
  
  others <- as.data.frame(findOverlaps(cranges, oranges, type = "within"))
  others <- others %>% 
    mutate(ctcf = chr_ctcfs$id[queryHits], sequence = chr_other$id[subjectHits], type = "other")
  
  hs <- bind_rows(promoters, genes, others) %>% select(ctcf, sequence, type)
  return(list(ctcfs = chr_ctcfs, hits = hs))
})

class_ctcfs <- ldply(lapply(ctcfs_hits, "[[", "ctcfs"))
ctcfs_hits <- ldply(lapply(ctcfs_hits, "[[", "hits"))

write_tsv(class_ctcfs, path = "data/class_ctcfs.tsv")
write_tsv(ctcfs_hits, path = "data/ctcfs_hits.tsv")

