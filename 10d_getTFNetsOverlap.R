library(igraph)
library(readr)
library(dplyr)
library(org.Hs.eg.db)

communities <- c(169,281,146,158,183,197,349,316,230,57,36,14,134,4,475)

tfs_files <- list.files("data/tfs/gtrd/promoters_bs/", pattern = "[a-zA-Z0-9]*.txt")
tfs_files <- tfs_files[sapply(paste0("data/tfs/gtrd/promoters_bs/", tfs_files), file.size) > 24]

interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
colnames(interactions)[1:2] <- c("from", "to")

net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)

all_coms <- lapply(communities, function(comm){
  v_in_comm <- vertices %>% filter(community == comm) %>% dplyr::select(ensemblID) %>% unlist(use.names = F)
  
  v_uniprots <- mapIds(org.Hs.eg.db,
                              keys = v_in_comm,
                              column="UNIPROT",
                              keytype="ENSEMBL",
                              multiVals="list")

  v_uniprots <- data.frame(ensembl_id = strtrim(names(unlist(v_uniprots)), 15), uniprot_id = unlist(v_uniprots))
  v_uniprots$file <- paste(v_uniprots$uniprot_id,".txt", sep = "")
  rownames(v_uniprots) <- NULL
  comm_files <- intersect(v_uniprots$file, tfs_files)
  v_uniprots %>% filter(file %in% comm_files) %>% write_tsv(paste0("data/tfs/gtrd/tfs_in_comm_", comm, ".txt"))
  comm_net <- induced.subgraph(net, v_in_comm)
  
  if(length(comm_files) > 0) {
    netss <- lapply(comm_files, function(cfile){
      tf_uniprot <- unlist(strsplit(cfile, "[.]"))[1]
      tf_ensembl <- v_uniprots %>% filter(file == cfile) %>% 
        dplyr::select(ensembl_id) %>% unlist()
      neigh <- read_tsv(paste0("data/tfs/gtrd/promoters_bs/", cfile), 
                        col_types = cols_only(to = col_character(), site_count = col_integer()),
                        col_names = c("to", "ensembl_id", "site_count"), skip = 1)
      if(nrow(neigh) > 0) {
        neigh$from <- tf_ensembl
        neigh <- neigh %>% dplyr::select(to, from, site_count)
        tf_net <- graph_from_data_frame(neigh, directed = F)
        
        return(get.data.frame(intersection(tf_net, comm_net), what = "edges"))
      }
    })
    netss <- bind_rows(netss)
    if(nrow(netss) > 0) {
      netss$community <- comm
      netss <- netss[!duplicated(netss %>% dplyr::select(to, from)), ]
      return(netss) 
    }
  }
  return(NULL)
})
bind_rows(all_coms) %>% write_tsv(paste0("data/tfs/gtrd/all_tfs_interactions.tsv"))

