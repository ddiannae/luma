library(igraph)
library(readr)
library(dplyr)
library(org.Hs.eg.db)

communities <- c(169,281,146,158,183,197,349,316,230,57,36,14,134,4,475)

tfs_files <- list.files("data/tfs/promoters_bs/", pattern = "[A-Z,0-9]*.txt")
biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt",   
                    col_names = c("ensemblID", "chr", "start", "end",  "gc", "type", "symbol"), skip = 1)
chrs <- c(as.character(1:22), "X")
biomart <- biomart %>% filter(chr %in% chrs)

interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
colnames(interactions)[1:2] <- c("from", "to")

net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)

all_coms <- lapply(communities, function(comm){
  v_in_comm <- vertices %>% filter(community == comm) %>% dplyr::select(ensemblID) %>% unlist()
  
  v_uniprots <- mapIds(org.Hs.eg.db,
                              keys = membership$ensemblID,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
  
  comm_files <- tfs_files[grep(pattern = comm, x = tfs_files)]
  comm_net <- induced.subgraph(net, v_in_comm)
  
  if(length(comm_files) > 0) {
    netss <- lapply(comm_files, function(cfile){
      neigh <- read_tsv(paste0("data/tfs/", cfile))
      tf <- strsplit(cfile, split = "_")[[1]][1]
      if(!tf %in% neigh$gene_symbol){
        neigh <- bind_rows(neigh, tibble(gene_symbol = tf, i_am_tf = 1, gene_description = tf, 
                                         chromosome_name = "",  upstream = 0, upstream_source = "",
                                         downstreamn =  nrow(neigh), organism = "Homo sapiens", p_value = 0))
      }
      neigh <- neigh %>% select(gene_symbol, i_am_tf) %>% 
        left_join(biomart, by = c("gene_symbol" = "symbol")) %>%
        select(ensemblID, everything()) %>% filter(!is.na(ensemblID))
      tf_int <- tibble(from = neigh %>% filter(gene_symbol == tf) %>% 
                         select(ensemblID) %>% unlist(use.names = F),
                       to = neigh$ensemblID)
      tf_net <- graph_from_data_frame(tf_int, directed = F, vertices = neigh)
      return(get.data.frame(intersection(tf_net, comm_net), what = "edges"))
    })
    return(netss)
  }
  return(NULL)
})

