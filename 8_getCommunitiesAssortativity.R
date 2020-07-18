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
    bynode <- rowSums(l_times_A) / rowSums(A)
    bynode_frac <- mean(bynode)
    total_frac <- sum(l_times_A)/sum(A)
    return(list(community_id = com, newman = newman_a, bynode = bynode_frac, totalfrac = total_frac))
  })
  comm_assort <- bind_rows(comm_assort)
  return(comm_assort)
}
 
n <- 20127  

#### luma
interactions <- read_tsv(file = paste0("data/network-tables/luma-", n, "-interactions.tsv"))
vertices <- read_tsv(file = paste0("data/network-tables/luma-", n, "-vertices.tsv"))
vertices <- vertices %>% mutate(coef = if_else(FDR > 0.05, 0, coef))

vertices <- vertices %>% mutate(exp = case_when(coef > 0 ~ "up", 
                                                coef < 0 ~ "down", 
                                                coef == 0 ~ "NA"))
vertices$exp <- as.factor(vertices$exp)

g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  
exp_assortativity <- getAssortativityByAttr(g, "exp")
exp_assortativity <- exp_assortativity %>% 
  inner_join(vertices %>% group_by(community) %>% summarise(mean_exp = mean(coef)), 
             by = c("community_id" = "community"))
write_tsv(exp_assortativity, path = "data/assortativity/luma-exp-assortativity.tsv")

chr_assortativity <- getAssortativityByAttr(g, "chr")
write_tsv(chr_assortativity, path = "data/assortativity/luma-chr-assortativity.tsv")

### healhty
interactions <- read_tsv(file = paste0("data/network-tables/healthy-", n, "-interactions.tsv"))
vertices <- read_tsv(file = paste0("data/network-tables/healthy-", n, "-vertices.tsv"))
g <- graph_from_data_frame(interactions, vertices = vertices, directed = FALSE)  
chr_assortativity <- getAssortativityByAttr(g, "chr")
write_tsv(chr_assortativity, path = "data/assortativity/healthy-chr-assortativity.tsv")
