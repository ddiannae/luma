library(readr)
library(tidyr)
library(dplyr)
library(igraph)
library(clusterProfiler)

algrthm = "multi_level"
conds <- c("healthy", "luma")

all_enr <-  lapply(conds, function(cond) {
  
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-20127-interactions.tsv"))
  vertices <- read_tsv(paste0("data/network-tables/", cond, "-20127-vertices.tsv"))
  colnames(interactions)[1:2] <- c("from", "to")
  universe <- vertices$ensemblID
  net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)
  gene_comm <- read_tsv(paste0("data/tfs/gtrd/", cond , "-tfs_interactions-", algrthm, ".tsv"))
  
  tf_inter <- lapply(unique(gene_comm$community), function(comm){
    comm_tf_inter <- read_tsv(paste0("data/tfs/gtrd/", cond, "-all_tfs_interactions-", comm, ".tsv"))
    
  })
  tf_inter <- bind_rows(tf_inter)
  tf_inter <- tf_inter %>% left_join(vertices %>% select(ensemblID, symbol), by=c("from"= "ensemblID")) %>%
    select(symbol, to, from)
  
  
  all_coms <- parallel::mclapply(X = unique(gene_comm$community), mc.cores = 5, FUN = function(comm){
  
    enrichs <- lapply(unique(tf_inter$from), function(tf){
      
      geneList <- unlist(neighbors(net, tf)$name)
      
      etf <- enricher(gene  = geneList,
                          universe  = universe,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2,
                          minGSSize     = 1,
                          maxGSSize = 10000,
                          TERM2GENE     = tf_inter)
      etf <- as.data.frame(etf)
      return(etf)
    })
      
    enrichs <- bind_rows(enrichs)
      
      if(nrow(enrichs) > 0) {
        enrichs$community <- comm
        return(enrichs) 
      }
  })
  all_coms <- bind_rows(all_coms) 
  all_coms$cond <- cond
  return(all_coms)
})
all_enr <- bind_rows(all_enr)
write_tsv(all_enr, file = "data/tfs/gtrd/all-tf-enrichments.tsv")
