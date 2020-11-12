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

healthy_interactions <- read_tsv("data/network-tables/healthy-20127-interactions.tsv")
healthy_interactions <- healthy_interactions %>% filter(interaction_type != "Trans")
healthy_vertices <- read_tsv("data/network-tables/healthy-vertices.tsv", 
                          col_types = cols(
                            chr = col_character()
                          ))
healthy_comm <- read_tsv(file = paste0("data/communities/healthy-communities.tsv"))

luma_net <- graph_from_data_frame(d = luma_interactions, vertices = luma_vertices, directed = FALSE)
healthy_net <- graph_from_data_frame(d = healthy_interactions, vertices = healthy_vertices, directed = FALSE)

getIntraCommunities <- function(net, membership) {
  
  intra_comm <- mclapply(X = unique(membership$community), 
           mc.cores = 75,
           FUN = function(n_comm){
             gene_list <- unlist(membership %>% filter(community == n_comm) %>%
                                  select(ensemblID))
               comm <- induced_subgraph(net, gene_list)
               if(is_connected(comm)){
                 return(n_comm)
               }  
             
             return(NULL)
           })
  return(na.omit(unlist(intra_comm)))
}

luma_intra_comm <- getIntraCommunities(luma_net, luma_comm)
healthy_intra_comm <- getIntraCommunities(healthy_net, healthy_comm)

getCommunitiesDiameter <- function(vertices, membership, communities) {
  diameters <- mclapply(X = unique(communities), 
                         mc.cores = 75,
                         FUN = function(n_comm){
                           gene_list <- unlist(membership %>% filter(community == n_comm) %>%
                                                 select(ensemblID))
                           comm_ver <- vertices %>% filter(ensemblID %in% gene_list)
                           return(list(community_id = n_comm, diameter = max(comm_ver$end) - min(comm_ver$start)))
                         })
  return(bind_rows(diameters))
} 

luma_intra_diam <- getCommunitiesDiameter(luma_vertices, luma_comm, luma_intra_comm)
healthy_intra_diam <- getCommunitiesDiameter(luma_vertices, luma_comm, healthy_intra_comm)

png(paste0("figures/communities/luma-intra-comm-diameter.png"), width = 1000, height = 500)
ggplot(luma_intra_diam, aes(x = diameter)) +
  geom_density() +
  scale_x_log10() +
  theme_bw()
dev.off()

png(paste0("figures/communities/healthy-intra-comm-diameter.png"), width = 1000, height = 500)
ggplot(healthy_intra_diam, aes(x = diameter)) +
  geom_density() +
  scale_x_log10() +
  theme_bw()
dev.off()

quantile(luma_intra_diam$diameter)
#0%         25%         50%         75%        100% 
#3491.0     60612.5    160807.0    521243.0 192474079.0 

quantile(healthy_intra_diam$diameter)
# 0%       25%       50%       75%      100% 
# 9343     51076    177629   4414307 220188132 

luma_intra_vertices <- luma_vertices %>% semi_join(luma_comm %>% 
                              filter(community %in% luma_intra_comm), by = "ensemblID") %>%
                              mutate(type = "gene")
healthy_intra_vertices <- healthy_vertices %>% semi_join(healthy_comm %>% 
                              filter(community %in% healthy_intra_comm), by = "ensemblID") %>%
                              mutate(type = "gene")

write_tsv(luma_intra_vertices, "data/communities/luma-intra-vertices.tsv")
write_tsv(luma_intra_diam, "data/communities/luma-intra-communities.tsv" )
write_tsv(healthy_intra_vertices, "data/communities/healthy-intra-vertices.tsv")
write_tsv(healthy_intra_diam, "data/communities/healthy-intra-communities.tsv" )
