library(readr)
library(dplyr)
library(parallel)
conds <- c("healthy", "luma")
algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector", "multi_level")

alg_pairs <- as.data.frame(combn(algorithms, 2))

m <- lapply(conds, function(cond){
  print(cond)
  
  m <- apply(alg_pairs, MARGIN = 2, FUN=function(ap)  {
    alg1 <- ap[1]
    alg2 <- ap[2]
    alg1_comms <- read_tsv(paste0("data/communities/", cond,
                               "-communities-", alg1,".tsv"))
    alg2_comms <-  read_tsv(paste0("data/communities/", cond,
                                "-communities-", alg2,".tsv")) 
  
  jaccs <- mclapply(unique(alg1_comms$community),
                    mc.cores = 5,
                    FUN = function(comm1) {
    genes1 <- alg1_comms %>% filter(community == comm1) %>%
      select(ensemblID) %>% unlist(use.names = F)
    
    jr <- lapply(unique(alg2_comms$community), function(comm2){
      genes2  <- alg2_comms %>% filter(community == comm2) %>%
        select(ensemblID) %>% unlist(use.names = F)
      return(list(comm1 = comm1, comm2 = comm2, genes1 = length(genes1), 
                  genes2 = length(genes2),
                  jaccard = length(intersect(genes1, genes2))/length(union(genes1, genes2))))
    })
    jr <- bind_rows(jr)
    return(jr)
  })
  
  jaccs <- bind_rows(jaccs)
  write_tsv(jaccs, paste0("data/communities/", cond, "-", alg1, "-", alg2, ".tsv"))
})
  
})