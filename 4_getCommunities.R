### Script to find communities (hot-spots) in the subnetworks spanned by 
### intra-chromosomal interactions. 
### Requires the type/intra/cond-all-distance-mi.txt file from the 
### 1_getIntraInteractions.R script
library(igraph)
library(readr)
library(dplyr)
library(ggplot2)

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

functions <-  c(cluster_leading_eigen, cluster_fast_greedy, cluster_louvain, 
                cluster_infomap)

m <- lapply(functions, function(f) {
  comms <- lapply(conds, function(cond) {
      
      interactions <- read_tsv(file = paste0("data/network-tables/", cond, "-", n, 
                                             "-interactions.tsv"))
    
      vertices <- read_tsv(file = paste0("data/network-tables/", cond, "-", n,
                                         "-vertices.tsv"))
      
      colnames(interactions)[1:3] <- c("from", "to", "weight")
      net <- graph_from_data_frame(interactions, directed=F, vertices = vertices)
      
      comm <- f(graph = net)
      algrthm <- gsub(" ", "_", comm$algorithm)
      print(paste0("Testing function: ", algrthm))
      names(comm$membership) <- comm$names
      
      df_comm <- data.frame(comm$names, comm$membership)
      colnames(df_comm) <- c("ensemblID", "community")
      
      vertices <- vertices %>% inner_join(df_comm, by = "ensemblID")
      

      write.table(df_comm, paste0("data/communities/", cond,  "-communities-", 
                                  algrthm, ".tsv"), quote = F, row.names = F,
                  col.names = T, sep = "\t")

      comm_info <- getComInfo(comm$membership, net)
      comm_info$cond <- cond 
      write.table(comm_info, paste0("data/communities/", cond, "-communities-info-",
                                    algrthm, ".tsv"), quote = F, row.names = F, 
                  col.names = T, sep = "\t")
      
      return(list(cond = cond, modularity = modularity(comm), algorithm = algrthm))
  })
  mod <- bind_rows(comms)
  return(mod)
})
m <- bind_rows(m)
m
# # A tibble: 8 x 3
# cond    modularity algorithm          
# <chr>        <dbl> <chr>              
# 1 healthy      0.696 leading_eigenvector
# 2 luma         0.892 leading_eigenvector
# 3 healthy      0.703 fast_greedy        
# 4 luma         0.934 fast_greedy        
# 5 healthy      0.752 multi_level        
# 6 luma         0.935 multi_level        
# 7 healthy      0.674 infomap            
# 8 luma         0.907 infomap  