library(readr)
library(dplyr)
library(plyr)
library(stringr)

annot <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                  col_names = c("ensembl_id", "chr", "start", "end", "karyoband"),
                  col_types = cols(col_character(), col_character(), col_integer(), 
                                   col_integer(), col_character()),
                  skip = 1)
genes <- read_tsv("data/genes_in_exp_matrix.txt", col_names = "ensembl_id")
chrs <- as.character(c(seq(1:22), "X"))

annot <- annot %>% semi_join(genes)
ctcfs <- read_tsv("data/GSE85108/GSM2257816_CTCF_E0h_MCF7_peaks.bed", 
                  col_names = c("chr", "start", "end", "id", "score"))
ctcfs <- ctcfs %>% mutate(chr = str_replace(chr, "chr", "")) %>% 
  filter(chr %in% chrs)
                            
nov_ctcfs <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  gs_in_c <- annot %>% filter(chr == ch)
  cs_in_c <- ctcfs %>% filter(chr == ch)
  
  cs <- lapply(1:nrow(cs_in_c), function(row){
    row <- cs_in_c[row, ]
    ends <- gs_in_c$end >= row["start"] & gs_in_c$end <= row["end"]
    starts <- gs_in_c$start >= row["start"] & gs_in_c$start <= row["end"]
    if(!(any(ends) | any(starts))) {
      return(row)
    }
    return(NULL)
  })
  cs <- ldply(compact(cs))
  return(cs)
})

nov_ctcfs <- ldply(nov_ctcfs)
                     
write_tsv(nov_ctcfs, path = "data/non_gene_overlaping_ctcfbs.tsv")
