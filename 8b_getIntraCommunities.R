library(readr)
library(dplyr)
library(igraph)
library(parallel)
library(ggplot2)

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

luma_interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
luma_interactions <- luma_interactions %>% filter(interaction_type != "Trans")
luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols(
                            chr = col_character()
                          ))


healthy_interactions <- read_tsv("data/network-tables/healthy-20127-interactions.tsv")
healthy_interactions <- healthy_interactions %>% filter(interaction_type != "Trans")
healthy_vertices <- read_tsv("data/network-tables/healthy-vertices.tsv", 
                          col_types = cols(
                            chr = col_character()
                          ))


algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector", "multi_level")

m <- lapply(algorithms, function(algrthm) {
  
  healthy_comm <- read_tsv(file = paste0("data/communities/healthy-communities-", algrthm, ".tsv"))
  luma_comm <- read_tsv(file = paste0("data/communities/luma-communities-", algrthm,".tsv"))
  
  luma_net <- graph_from_data_frame(d = luma_interactions, vertices = luma_vertices, directed = FALSE)
  healthy_net <- graph_from_data_frame(d = healthy_interactions, vertices = healthy_vertices, directed = FALSE)
  
  healthy_chr_assort <- read_tsv(file = paste0("data/assortativity/healthy-chr-assortativity-", algrthm, ".tsv"))
  luma_chr_assort <- read_tsv(file = paste0("data/assortativity/luma-chr-assortativity-", algrthm,".tsv"))
  
  luma_intra_comm <- luma_chr_assort %>% filter(diffraction == 1) %>% 
    select(community_id) %>% unlist(use.names = F)
    
  healthy_intra_comm <- healthy_chr_assort%>% filter(diffraction == 1) %>% 
    select(community_id) %>% unlist(use.names = F)
  
  luma_intra_diam <- getCommunitiesDiameter(luma_vertices, luma_comm, luma_intra_comm)
  healthy_intra_diam <- getCommunitiesDiameter(luma_vertices, luma_comm, healthy_intra_comm)
  
  luma_intra_vertices <- luma_vertices %>% semi_join(luma_comm %>% 
                                                       filter(community %in% luma_intra_comm), by = "ensemblID") %>%
    mutate(type = "gene")
  healthy_intra_vertices <- healthy_vertices %>% semi_join(healthy_comm %>% 
                                                             filter(community %in% healthy_intra_comm), by = "ensemblID") %>%
    mutate(type = "gene")
  
  write_tsv(luma_intra_vertices, paste0("data/communities/luma-intra-vertices-", algrthm, ".tsv"))
  write_tsv(luma_intra_diam, paste0("data/communities/luma-intra-communities-", algrthm, ".tsv" ))
  write_tsv(healthy_intra_vertices, paste0("data/communities/healthy-intra-vertices-", algrthm, ".tsv"))
  write_tsv(healthy_intra_diam, paste0("data/communities/healthy-intra-communities-", algrthm, ".tsv" ))
  
})

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

conds <- c("healthy", "luma")
all_comms <- lapply(conds, function(cond) {
  m <- lapply(algorithms, function(algrthm) {
    comms <- read_tsv(file = paste0("data/communities/", cond, "-communities-info-", algrthm, ".tsv"),
                      col_types = cols(chr = col_character() ))
    intra_comms <- read_tsv(file = paste0("data/communities/",cond, "-intra-communities-", algrthm, ".tsv"))
    comms$algrthm <- algrthm
    comms$cond <- cond
    comms <- comms %>% mutate(community_type = 
                                            ifelse(com_id %in% intra_comms$community_id, 
                                                   "intra", "inter"))
    return(comms)
    
    })
  m <- bind_rows(m)
})

all_comms <- bind_rows(all_comms)
all_comms %>% group_by(cond, algrthm, community_type) %>% tally()

# # A tibble: 16 x 4
# # Groups:   cond, algrthm [8]
# cond    algrthm             community_type     n
# <chr>   <chr>               <chr>          <int>
# 1 healthy fast_greedy         inter            325
# 2 healthy fast_greedy         intra             75
# 3 healthy infomap             inter            768
# 4 healthy infomap             intra             83
# 5 healthy leading_eigenvector inter            283
# 6 healthy leading_eigenvector intra             71
# 7 healthy multi_level         inter            291
# 8 healthy multi_level         intra             71
# 9 luma    fast_greedy         inter             87
# 10 luma    fast_greedy         intra            614
# 11 luma    infomap             inter             93
# 12 luma    infomap             intra            826
# 13 luma    leading_eigenvector inter             84
# 14 luma    leading_eigenvector intra            594
# 15 luma    multi_level         inter             87
# 16 luma    multi_level         intra            614

all_comms %>% filter(order >= 5) %>% group_by(cond, algrthm, community_type) %>% tally()
# cond    algrthm             community_type     n
# <chr>   <chr>               <chr>          <int>
# 1 healthy fast_greedy         inter             50
# 2 healthy infomap             inter            386
# 3 healthy infomap             intra              1
# 4 healthy leading_eigenvector inter             32
# 5 healthy leading_eigenvector intra              1
# 6 healthy multi_level         inter             41
# 7 luma    fast_greedy         inter             40
# 8 luma    fast_greedy         intra             77
# 9 luma    infomap             inter             39
# 10 luma    infomap             intra            194
# 11 luma    leading_eigenvector inter             37
# 12 luma    leading_eigenvector intra             58
# 13 luma    multi_level         inter             40
# 14 luma    multi_level         intra             77