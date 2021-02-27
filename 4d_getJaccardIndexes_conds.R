library(readr)
library(dplyr)
library(parallel)

conds <- c("healthy", "luma")
algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector", "multi_level")

m <- lapply(algorithms, function(alg) {
 
    healthy_comms <- read_tsv(paste0("data/communities/", conds[1],
                                  "-communities-", alg,".tsv"))
    luma_comms <-  read_tsv(paste0("data/communities/", conds[2],
                                   "-communities-", alg,".tsv")) 
    
    jaccs <- mclapply(unique(healthy_comms$community),
                      mc.cores = 5,
                      FUN = function(comm1) {
                        genes1 <- healthy_comms %>% filter(community == comm1) %>%
                          select(ensemblID) %>% unlist(use.names = F)
                        
                        jr <- lapply(unique(luma_comms$community), function(comm2){
                          genes2  <- luma_comms %>% filter(community == comm2) %>%
                            select(ensemblID) %>% unlist(use.names = F)
                          
                          return(list(healthy = comm1, luma = comm2, healthy_n = length(genes1), 
                                      luma_n = length(genes2),
                                      jaccard = length(intersect(genes1, genes2))/length(union(genes1, genes2))))
                        })
                        jr <- bind_rows(jr)
                        return(jr)
                      })
    
    jaccs <- bind_rows(jaccs)
    write_tsv(jaccs, paste0("data/communities/healthy-luma-", alg, ".tsv"))

})