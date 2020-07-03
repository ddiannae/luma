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

ctcfs <- read_tsv("data/GSE85108/GSM2257816_CTCF_E0h_MCF7_peaks_hg38.bed", 
                  col_names = c("chr", "start", "end", "id", "score"))
ctcfs <- ctcfs %>% mutate(chr = str_remove(chr, "chr")) %>% 
  filter(chr %in% chrs)

promoters <- protein_coding %>% mutate(end = start + 5000, start = start - 5000)
genes <- protein_coding %>% mutate(start = start + 5000)
genes <- genes %>% filter(end > start)

class_ctcfs <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
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
                            ohits = countOverlaps(cranges, oranges, type = "within"),
                            promoter = unlist(chr_promoters[findOverlaps(cranges, pranges, type = "within", select = "first"), "id"]),
                            gen = unlist(chr_genes[findOverlaps(cranges, granges, type = "within", select = "first"), "id"]),
                            other = unlist(chr_other[findOverlaps(cranges, oranges, type = "within", select = "first"), "id"]))
  return(chr_ctcfs)
})

  
class_ctcfs <- ldply(class_ctcfs)
write_tsv(class_ctcfs, path = "data/all_ctcfs.tsv")
