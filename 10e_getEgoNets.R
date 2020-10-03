library(igraph)
library(readr)

communities <- c(169,281,146,158,183,197,349,316,230,57,36,14,134,4,475)
interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
colnames(interactions)[1:2] <- c("from", "to")
net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)

tf_interactions <- read_tsv("data/tfs/gtrd/all_tfs_interactions.tsv")

igraph_options(vertex.size = 5, vertex.label.color = "black", 
               vertex.frame.color = "azure", vertex.color = "wheat2", vertex.alpha =0.5,
               edge.color = "gray90")

all_coms <- lapply(communities, function(comm){
  v_in_comm <- vertices %>% filter(community == comm) %>% 
    dplyr::select(ensemblID) %>% unlist(use.names = F)
  
  comm_tfs <- read_tsv(paste0("data/tfs/gtrd/tfs_in_comm_", comm, ".txt"))
  comm_net <- induced.subgraph(net, v_in_comm)
  comm_inter <- tf_interactions %>% filter(community == comm)
  
  gs <- lapply(comm_tfs$ensembl_id, function(tf) {
   
    eg <-  make_ego_graph(comm_net,
                          order = 1,
                          nodes = tf)[[1]]
    eg_inter <- comm_inter %>% filter(to %in% V(eg)$name & from %in% V(eg)$name)
    veg_inter <- unique(c(tf, unlist(eg_inter$from), unlist(eg_inter$to)))

    ### Edge colors
    # All edges
    E(eg)$color <- "grey"
    # Regulatory edges
    if(nrow(eg_inter) > 0) {
      for(i in seq_along(eg_inter$from)) {
        E(eg)[eg_inter$from[i] %--% eg_inter$to[i]]$color <- "tomato4"
      }
    }
   
    ## Node labels
    labels <- ifelse(V(eg)$name %in% veg_inter | V(eg)$name %in% comm_tfs$ensembl_id, V(eg)$symbol, NA)
    ## Node colors
    colors <- ifelse(V(eg)$name %in% comm_tfs$ensembl_id,  "paleturquoise3", "wheat2")
    ## Node Sizes
    sizes <- ifelse(V(eg)$name %in% comm_tfs$ensembl_id,  15, 5)

    png(paste0("figures/tfs/", tf, ".png"))
      plot(eg, vertex.color = colors, layout = layout_with_fr,
         vertex.label = labels,
         vertex.label.family="Helvetica", vertex.size = sizes
         )
    dev.off()
  })
})
  
