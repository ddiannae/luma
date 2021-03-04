library(readr)
library(janitor)
library(dplyr)

algorithms <- c("fast_greedy", "infomap" ,"leading_eigenvector", "multi_level")

all_enrichments <- lapply(algorithms, function(alg){
    enrich_sanos <- read_tsv(paste0("data/enrich/healthy-go-enrichments-", 
                                    alg,".tsv")) %>%
                  clean_names()
    enrich_luma <-  read_tsv(paste0("data/enrich/luma-go-enrichments-", 
                                    alg,".tsv")) %>%
                    clean_names()
    
    jaccs <- lapply(unique(enrich_sanos$commun), function(comm1) {
      terms1 <- enrich_sanos %>% filter(commun == comm1) %>%
        select(id) %>% unlist(use.names = F)
      
      jr <- lapply(unique(enrich_luma$commun), function(comm2){
        terms2  <- enrich_luma %>% filter(commun == comm2) %>%
          select(id) %>% unlist(use.names = F)
        return(list(comm1 = comm1, comm2 = comm2, terms1 = length(terms1), 
                    terms2 = length(terms2),
                    jaccard = length(intersect(terms1, terms2))/length(union(terms1, terms2))))
      })
      jr <- bind_rows(jr)
      return(jr)
    })
    
    jaccs <- bind_rows(jaccs)
    write_tsv(jaccs, paste0("data/enrich/luma-healthy-", alg, ".tsv"))
})
