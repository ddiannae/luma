library(readr)
library(dplyr)

n <- 20127  
conds <- c("healthy", "luma")

### First save the vertices with communities.
m <- lapply(conds, function(cond) {
  interactions <- read_tsv(file = paste0("data/network-tables/", cond, "-", n, "-interactions.tsv"))
  
  vertices <- read_tsv(file = paste0("data/network-tables/", cond, "-", n, "-vertices.tsv"))
  
  ### Adding community info
  gene_sets <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  vertices <- vertices %>% inner_join(gene_sets, by = "ensemblID")
  
  ### Adding DEG info
  if(cond != "healthy") {
    deg <- read_tsv(file = paste0("data/", cond,  "-deg-ebayes.tsv"))
    vertices <-  vertices %>% inner_join(deg, by = c("ensemblID" = "gene"))
  }
  
  write.table(vertices, file = paste0("data/network-tables/", cond, "-", n, "-vertices.tsv") , 
              quote = F, row.names = F, col.names = T, sep = "\t")
  
})
    