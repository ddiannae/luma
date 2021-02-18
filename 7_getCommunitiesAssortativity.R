library(readr)
library(dplyr)
library(igraph)

getAssortativityByAttr <- function(gr, n_attr) {
  vs <- as_data_frame(gr, what = "vertices")
  comm_assort <- lapply(unique(vs$community), function(com){
    vc <- vs %>% filter(community == com) %>% 
      select(name) %>% unlist()
    com_graph <- induced_subgraph(gr, vc)
    v_attr <- as_data_frame(com_graph, what = "vertices")[, n_attr]
    newman_a <- assortativity_nominal(com_graph,
                                      types = as.numeric(as.factor(v_attr)), 
                                      directed = T)
    A <- as.matrix(get.adjacency(com_graph))
    l_times_A <- outer(v_attr, v_attr, `==`) * A
    l_dif_A <- outer(v_attr, v_attr, `!=`) * A
    bynode <- rowSums(l_times_A) / rowSums(A)
    bynode_frac <- mean(bynode)
    total_frac <- sum(l_times_A)/sum(A)
    dif_fraction <- (sum(l_times_A)-sum(l_dif_A))/sum(A)
    return(list(community_id = com, newman = newman_a, bynode = bynode_frac, 
                totalfrac = total_frac, diffraction=dif_fraction))
  })
  comm_assort <- bind_rows(comm_assort)
  return(comm_assort)
}
 
n <- 20127  
algorithms <- c("leading_eigenvector", "fast_greedy", "multi_level", "infomap")

### Luma
interactions <- read_tsv(file = paste0("data/network-tables/luma-", n, "-interactions.tsv"))
vertices <- read_tsv(file = paste0("data/network-tables/luma-", n, "-vertices.tsv"))
vertices <- vertices %>% select(-community)
# vertices <- vertices %>% mutate(exp = case_when(lfc > 0 ~ "up", 
#                                                 lfc < 0 ~ "down", 
#                                                 lfc == 0 ~ "NA"))
# vertices$exp <- as.factor(vertices$exp)

# g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  
# exp_assortativity <- getAssortativityByAttr(g, "exp")
# exp_assortativity <- exp_assortativity %>% 
#   inner_join(vertices %>% group_by(community) %>% 
#                summarise(mean_diff_exp = mean(lfc),
#                          mean_avg_exp  =  mean(avg_exp),
#                          mean_avg_log2_exp = mean(avg_log2_exp),
#                          m <- mapply(function(cond, algrthm) {
#                            
#                            cat("Working with condition: ", cond, ", ", algrthm, "\n")
#                            
#                            membership <- read_tsv(file = paste0("data/communities/", cond,  "-communities-",
#                                                                 algrthm, ".tsv"))
#                            mean_avg_cpm_exp = mean(avg_cpm_exp), 
#                          mean_avg_log2_cpm_exp = mean(avg_log2_cpm_exp)),
#              by = c("community_id" = "community"))
# write_tsv(exp_assortativity, path = "data/assortativity/luma-exp-assortativity.tsv")

m <- lapply(algorithms, function(algrthm) {
  
  cat("Working with algorithm:", algrthm, "\n")
  
  membership <- read_tsv(file = paste0("data/communities/luma-communities-",
                                       algrthm, ".tsv"))
  vertices <- vertices %>% left_join(membership)
  g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  
  chr_assortativity <- getAssortativityByAttr(g, "chr")
  write_tsv(chr_assortativity, path = paste0("data/assortativity/luma-chr-assortativity-",algrthm,".tsv"))
})

## Healthy
interactions <- read_tsv(file = paste0("data/network-tables/healthy-", n, "-interactions.tsv"))
vertices <- read_tsv(file = paste0("data/network-tables/healthy-", n, "-vertices.tsv"))
vertices <- vertices %>% select(-community)

m <- lapply(algorithms, function(algrthm) {
  
  cat("Working with algorithm:", algrthm, "\n")
  
  membership <- read_tsv(file = paste0("data/communities/healthy-communities-",
                                       algrthm, ".tsv"))
  vertices <- vertices %>% left_join(membership)
  g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  
  chr_assortativity <- getAssortativityByAttr(g, "chr")
  write_tsv(chr_assortativity, path = paste0("data/assortativity/healthy-chr-assortativity-",algrthm,".tsv"))
})
