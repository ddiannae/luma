library(readr)
library(dplyr)
library(igraph)

n <- 20127  
interactions <- read_tsv(file = paste0("data/network-tables/luma-", n, "-interactions.tsv"))
vertices <- read_tsv(file = paste0("data/network-tables/luma-", n, "-vertices.tsv"))

g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  

lapply(unique(vertices$community), function(com){
  vc <- vertices %>% filter(community == com) %>% select(ensemblID) %>% unlist()
  community <- induced_subgraph(g, vc)
  
})