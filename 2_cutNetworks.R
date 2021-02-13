library(readr)
library(dplyr)

conds <- c("healthy", "luma")

all_interactions <- lapply(conds, function(cond){
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-interactions.tsv"))
  interactions$cond <- cond
  return(interactions)
})

all_interactions <- bind_rows(all_interactions)
all_interactions %>% group_by(cond) %>% tally()

#1 healthy 20127
#2 luma    30269

n <- all_interactions %>% group_by(cond) %>% tally() %>% 
  filter(n == min(n)) %>% select(n) %>% as.integer()

lapply(conds, function(cond){
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-interactions.tsv"))
  interactions <- interactions %>% arrange(desc(MI)) 
  interactions <- head(interactions, n = n)
  interactions$row_num <- 1:nrow(interactions)
  
  vertices <- read_tsv(paste0("data/network-tables/", cond, "-vertices.tsv"))
  vertices <- vertices %>% filter(ensemblID %in% unique(c(interactions$source, interactions$target)))
  
  write.table(vertices, file = paste0("data/network-tables/", cond, "-", n, "-vertices.tsv") , 
              quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(interactions, file = paste0("data/network-tables/", cond, "-", n, "-interactions.tsv"),
              quote = F, row.names = F, col.names = T, sep = "\t")
})
