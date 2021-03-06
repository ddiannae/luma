library(readr)
library(tidyr)
library(dplyr)

comms <- c(208, 54, 572)
interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/communities/luma-communities-multi_level.tsv")
tf_interactions <- read_tsv("data/tfs/gtrd/all_tfs_interactions-multi_level.tsv")

m <- lapply(comms, function(comm){
  v_in_comm <- vertices %>% filter(community == comm)
  comm_interactions <- interactions %>%
    filter(source %in% v_in_comm$ensemblID & target %in% v_in_comm$ensemblID)
  
  comm_interactions$isRegulatory = ifelse(comm_interactions$row_num %in% tf_interactions$row_num, 1, 0)
  write_tsv(comm_interactions, paste0("data/tfs/interactions-comm-", comm, ".tsv"))
})
