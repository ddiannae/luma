library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

annot <- read_tsv("../data/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                  col_names = c("ensemblID", "chr", "start", "end", "karyoband"),
                  col_types = cols(col_character(), col_character(), col_integer(), 
                                   col_integer(), col_character()),
                  skip = 1)
genes <- read_tsv("../data/genes_in_exp_matrix.txt", col_names = "ensemblID")

annot <- annot %>% semi_join(genes)
ctcfs <- read_csv("../data/ctcfbs.csv", 
                  col_names = c("ID", "specie", "chr", "cell", "exp", "source"), skip = 1)
ctcfs <- ctcfs %>% select("chr")
ctcfs <- ctcfs %>% separate(col = "chr", into = c("chr", "pos"), sep = ":")
ctcfs <- ctcfs %>% separate(col = "pos", into = c("start", "end"), sep = "-", convert = T)
ctcfs$chr <- str_replace(ctcfs$chr, "chr", "")

chrs <- as.character(c(seq(1:22), "X"))

nov_ctcfs <- parallel::mclapply(X = chrs, mc.cores = 5, FUN = function(c){
  gs_in_c <- annot %>% filter(chr == c)
  cs_in_c <- ctcfs %>% filter(chr == c)
  
  cs <- apply(cs_in_c, MARGIN = 1, function(row){
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
                     
write_tsv(nov_ctcfs, path = "../data/non_gene_overlaping_ctcfbs.tsv")
