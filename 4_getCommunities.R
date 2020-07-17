### Script to find communities (hot-spots) in the subnetworks spanned by 
### intra-chromosomal interactions. 
### Requires the type/intra/cond-all-distance-mi.txt file from the 
### 1_getIntraInteractions.R script
library(igraph)
library(readr)
library(dplyr)

getComInfo <- function(cmembership, network){
  comp.info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    com <- induced.subgraph(network, mem)
    prs <- page.rank(com)
    chrs <- table(V(com)$chr)
    return(data.frame(com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], order = length(V(com)), 
                      size = length(E(com))))
  })
  comp.info <- plyr::compact(comp.info)
  comp.info <- bind_rows(comp.info)
  comp.info <- comp.info %>% arrange(desc(size))
  return(comp.info)
}

n <- 20127  
conds <- c("healthy", "luma")

m <- lapply(conds, function(cond) {
  interactions <- read_tsv(file = paste0("data/network-tables/", cond, "-", n, "-interactions.tsv"))

  vertices <- read_tsv(file = paste0("data/network-tables/", cond, "-", n, "-vertices.tsv"))
  
  colnames(interactions)[1:2] <- c("from", "to")
  net <- graph_from_data_frame(interactions, 
  directed=F, vertices = vertiices)
  comm <- cluster_louvain(graph = net)
  names(comm$membership) <- comm$names
  
  df_comm <- data.frame(comm$names, comm$membership)
  colnames(df_comm) <- c("ensemblID", "community")
  
  vertices <- vertices %>% inner_join(df_comm, by = "ensemblID")
  
  write.table(vertices, file = paste0("data/network-tables/", cond, "-", n, "-vertices.tsv") , 
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  write.table(df_comm, paste0("data/communities/", cond,  "-communities.tsv"), 
            quote = F, row.names = F, col.names = T, sep = "\t")
  
  comm_info <- getComInfo(comm$membership, net)
  write.table(comm_info, paste0("data/communities/", cond, "-communities-info.tsv"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
})
