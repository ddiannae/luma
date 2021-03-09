library(igraph)
library(readr)
library(dplyr)
library(org.Hs.eg.db)

tfs_files <- list.files("data/tfs/gtrd/promoters_bs/", pattern = "[a-zA-Z0-9]*.txt")
tfs_files <- tfs_files[sapply(paste0("data/tfs/gtrd/promoters_bs/", tfs_files), file.size) > 24]


algrthm = "multi_level"
conds <- c("healthy", "luma")
all_enrh <-  lapply(conds, function(cond) {
  
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-20127-interactions.tsv"))
  vertices <- read_tsv(paste0("data/network-tables/", cond, "-20127-vertices.tsv"))
  colnames(interactions)[1:2] <- c("from", "to")
  
  net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)
  
  gene_comm <- read_tsv(paste0("data/communities/", cond , "-communities-", algrthm, ".tsv"))
    
  all_coms <- lapply(unique(gene_comm$community), function(comm){
    
    v_in_comm <- gene_comm %>% dplyr::filter(community == comm) %>% 
      dplyr::select(ensemblID) %>% unlist(use.names = F)
    
    v_uniprots <- mapIds(org.Hs.eg.db,
                                keys = v_in_comm,
                                column="UNIPROT",
                                keytype="ENSEMBL",
                                multiVals="list")
  
    v_uniprots <- data.frame(ensembl_id = strtrim(names(unlist(v_uniprots)), 15), uniprot_id = unlist(v_uniprots))
    v_uniprots$file <- paste(v_uniprots$uniprot_id,".txt", sep = "")
    rownames(v_uniprots) <- NULL
    comm_files <- intersect(v_uniprots$file, tfs_files)
    v_uniprots %>% filter(file %in% comm_files) %>% 
      write_tsv(paste0("data/tfs/gtrd/", cond, "-tfs_in_comm_", comm, "-", algrthm, ".txt"))
    comm_net <- induced.subgraph(net, v_in_comm)
    
    if(length(comm_files) > 0) {
      netss <- lapply(comm_files, function(cfile) {
        tf_uniprot <- unlist(strsplit(cfile, "[.]"))[1]
        
        tf_ensembl <- v_uniprots %>% filter(file == cfile) %>% 
          dplyr::select(ensembl_id) %>% unlist()
        
        neigh <- read_tsv(paste0("data/tfs/gtrd/promoters_bs/", cfile), 
                          col_types = cols_only(to = col_character(), site_count = col_integer()),
                          col_names = c("to", "ensembl_id", "site_count"), skip = 1)
        
        if(nrow(neigh) > 0) {
          neigh$from <- tf_ensembl
          neigh <- neigh %>% dplyr::select(from, to, site_count)
          tf_net <- graph_from_data_frame(neigh, directed = F)
          tf_inter <- igraph::as_data_frame(intersection(tf_net, comm_net), what = "edges")
          
          inter_not_tf <- sum(!neighbors(net, tf_ensembl)$name %in% neigh$to)
          tf_in_net <- length(intersect(neigh$to, vertices$ensemblID))
          
          if(nrow(tf_inter) > 0){
            tf_inter$tf = tf_ensembl
          }
          
          
          d <- data.frame(not_tf_neigh = c(tf_in_net - nrow(tf_inter), 
                                           nrow(vertices) - tf_in_net - inter_not_tf), 
                          tf_neigh = c(nrow(tf_inter), inter_not_tf), 
                          tf = c(tf_ensembl, tf_ensembl),
                          type = c("regulatory", "not_regulatory"))
     
          return(list(inter = tf_inter, hyper_df = d, all_tf_inter = neigh))
        }
      })
      
      reg_nets <- bind_rows(lapply(netss, "[[","inter" ))
      
      bind_rows(lapply(netss, "[[","hyper_df" )) %>%
        write_tsv(paste0("data/tfs/gtrd/", cond, "-fisher_matrices-", comm, ".tsv"))
      
      bind_rows(lapply(netss, "[[","all_tf_inter" )) %>% 
        dplyr::select(from, to) %>% 
        write_tsv(paste0("data/tfs/gtrd/", cond, "-all_tfs_interactions-", comm, ".tsv"))
      
      if(nrow(reg_nets) > 0) {
        reg_nets$community <- comm
        reg_nets <- reg_nets[!duplicated(reg_nets %>% dplyr::select(to, from, tf)), ]
        return(reg_nets) 
      }
    }
    return(NULL)
  })
  bind_rows(all_coms) %>% write_tsv(paste0("data/tfs/gtrd/", cond,"-tfs_interactions-", algrthm, ".tsv"))
})


##                 gene.not.interest gene.in.interest
## In_category                  //anotados como regulados en la red, que no son vecinos  76 // Vecinos de TF que están anotados como regulatorios
## not_in_category             // todos los demás en la red            // Vecinos de TF que no están anotados como regulatorios
## ----------------------------------------------------------------------------------------------------
## Red                                #vecinos de TF

