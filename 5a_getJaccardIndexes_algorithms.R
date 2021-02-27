library(readr)
library(janitor)
library(dplyr)

conds <- c("healthy", "luma")
algorithms <- c("fast_greedy", "infomap" ,"leading_eigenvector", "multi_level")

alg_pairs <- as.data.frame(combn(algorithms, 2)) %>%
  bind_cols(t(data.frame(algorithms, algorithms)))

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
        return(list(comm1 = comm1, comm2 = comm2, terms1 = length(terms1), 
                    terms2 = length(terms2),
                    jaccard = length(intersect(terms1, terms2))/length(union(terms1, terms2))))
      })
      jr <- bind_rows(jr)
      return(jr)
    })

    jaccs <- bind_rows(jaccs)
    write_tsv(jaccs, paste0("data/enrich/", cond, "-", alg1, "-", alg2, ".tsv"))
  })
})


all_enrichments <- lapply(conds, function(cond){
  algo_info <- lapply(algorithms, function(alg) {
    comm_info <- read_tsv(paste0("data/communities/", cond, "-communities-info-", 
                                 alg, ".tsv"))
    enrich <- read_tsv(paste0("data/enrich/", cond,
                               "-go-enrichments-", alg,".tsv")) %>%
      clean_names()
    enrich_comm <- enrich %>% group_by(commun) %>% tally()
    ncomm_plus5 <- comm_info %>% filter(order >= 5) %>% nrow()
    enrich_plus10 <- enrich_comm %>% filter(n >= 10) %>% nrow()
    return(list(algorithm = alg, n_comm = nrow(comm_info), n_comm_5 = ncomm_plus5,
                n_comm_enrich = nrow(enrich_comm), n_comm_enrich_10 = enrich_plus10,
                n_terms_enrich = nrow(enrich)))
  })
  algo_info <- bind_rows(algo_info)
  algo_info$cond <- cond
  return(algo_info)
})
all_info <- bind_rows(all_enrichments)
write_tsv(all_info, paste0("data/enrich/all-algorithms-comparison.tsv"))
