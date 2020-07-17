library(readr)
library(dplyr)
library(igraph)
library(parallel)
library(ggplot2)

luma_interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
luma_interactions <- luma_interactions %>% filter(interaction_type != "Trans")
luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols(
                            chr = col_character()
                          ))

luma_comm <- read_tsv(file = paste0("data/communities/luma-communities.tsv"))

luma_net <- graph_from_data_frame(d = luma_interactions, vertices = luma_vertices, directed = FALSE)

getIntraCommunities <- function(net, membership) {
  
  intra_comm <- mclapply(X = unique(membership$community), 
           mc.cores = 75,
           FUN = function(n_comm){
             gene_list <- unlist(membership %>% filter(community == n_comm) %>%
                                  select(ensemblID))
               comm <- induced_subgraph(luma_net, gene_list)
               if(is_connected(comm)){
                 return(n_comm)
               }  
             
             return(NULL)
           })
  return(na.omit(unlist(intra_comm)))
}

intra_comm <- getIntraCommunities(luma_net, luma_comm)

getCommunitiesDiameter <- function(vertices, membership, communities) {
  diameters <- mclapply(X = unique(communities), 
                         mc.cores = 75,
                         FUN = function(n_comm){
                           gene_list <- unlist(membership %>% filter(community == n_comm) %>%
                                                 select(ensemblID))
                           comm_ver <- vertices %>% filter(ensemblID %in% gene_list)
                           return(list(comm = n_comm, diam = max(comm_ver$end) - min(comm_ver$start)))
                         })
  return(bind_rows(diameters))
} 

intra_diam <- getCommunitiesDiameter(luma_vertices, luma_comm, intra_comm)

png(paste0("figures/communities/intra-comm-diameter.png"), width = 1000, height = 500)
ggplot(intra_diam, aes(x = diam)) +
  geom_density() +
  scale_x_log10() +
  theme_bw()
dev.off()

quantile(intra_diam$diam)
#0%         25%         50%         75%        100% 
#3491.0     60612.5    160807.0    521243.0 192474079.0 

intra_vertices <- luma_vertices %>% semi_join(luma_comm %>% 
                              filter(community %in% intra_comm), by = "ensemblID") %>%
                              mutate(type = "gene")

write_tsv(intra_vertices, path = "data/communities/luma-intra-vertices.tsv")
