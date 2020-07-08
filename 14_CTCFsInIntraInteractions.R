library(readr)
library(dplyr)
library(igraph)
library(parallel)
library(ggplot2)

luma_interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
luma_interactions <- luma_interactions %>% filter(interaction_type != "Trans")
luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(
                            ensemblID = col_character(),
                            chr = col_character(),
                            start = col_integer(),
                            end = col_integer()
                          ))

luma_comm <- read_tsv(file = paste0("data/luma-communities.tsv"))

luma_net <- graph_from_data_frame(d = luma_interactions, vertices = luma_vertices, directed = FALSE)

getIntraCommunities <- function(net, membership) {
  
  intra_comm <- mclapply(X = unique(membership$community), 
           mc.cores = 75,
           FUN = function(n_comm){
             gene_list <- unlist(membership %>% filter(community == n_comm) %>%
                                  select(ensemblID))
             #if(length(gene_list) > 5) {
               comm <- induced_subgraph(luma_net, gene_list)
               if(is_connected(comm)){
                 return(n_comm)
               }  
             #}
             
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
ggplot(intra_diam, aes(x = diam)) +
  geom_density() +
  scale_x_continuous(trans='log2') +
  theme_bw()

quantile(intra_diam$diam)
# 0%          25%          50%          75%         100% 
# 3491.00     72223.75    209088.50   2525232.75 246350418.00 
