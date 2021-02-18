library(readr)
library(janitor)
library(dplyr)

conds <- c("healthy", "luma")
algorithms <- c("leading_eigenvector", "fast_greedy", "multi_level", "infomap")

alg_pairs <- as.data.frame(combn(algorithms, 2))
all_enrichments <- lapply(conds, function(cond){
  m <- apply(alg_pairs, MARGIN = 2, FUN=function(ap)  {
    alg1 <- ap[1]
    alg2 <- ap[2]
    enrich1 <- read_tsv(paste0("data/enrich/", cond,
                                   "-go-enrichments-", alg1,".tsv")) %>%
              clean_names()
    enrich2 <-  read_tsv(paste0("data/enrich/", cond,
                                "-go-enrichments-", alg2,".tsv")) %>%
            clean_names()

    jaccs <- lapply(unique(enrich1$commun), function(comm1) {
      terms1 <- enrich1 %>% filter(commun == comm1) %>%
        select(id) %>% unlist(use.names = F)

      jr <- lapply(unique(enrich2$commun), function(comm2){
        terms2  <- enrich2 %>% filter(commun == comm2) %>%
          select(id) %>% unlist(use.names = F)
        return(list(comm1 = comm1, comm2 = comm2, jaccard =
                      length(intersect(terms1, terms2))/length(union(terms1, terms2))))
      })
      jr <- bind_rows(jr)
      return(jr)
    })

    jaccs <- bind_rows(jaccs)
    write_tsv(jaccs, paste0("data/enrich/", cond, "-", alg1, "-", alg2, ".tsv"))
  })
})
