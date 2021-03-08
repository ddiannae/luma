library(igraph)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)

interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
colnames(interactions)[1] = "from"
colnames(interactions)[2] = "to"
biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt",   
                    col_names = c("ensemblID", "chr", "start", "end",  "gc", "type", "symbol"), skip = 1)

algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector","multi_level")

m <- lapply(algorithms, function(algrthm) {
  
  communities <- read_tsv(paste0("data/enrich/communities-more-20-terms-", algrthm, ".tsv")) %>% unlist()
  gene_comm <- read_tsv(paste0("data/communities/luma-communities-", algrthm, ".tsv"))
  
  tf_interactions <- read_tsv(paste0("data/tfs/gtrd/all_tfs_interactions-", algrthm, ".tsv"))
  
  all_tfs <- lapply(communities, function(comm){
    v_in_comm <- gene_comm %>% filter(community == comm) %>% 
      dplyr::select(ensemblID) %>% unlist(use.names = F)
    
    comm_tfs <- read_tsv(paste0("data/tfs/gtrd/tfs_in_comm_", comm, "-", algrthm, ".txt"))
    e_in_comm <- interactions %>% filter(from %in% v_in_comm & to %in% v_in_comm)
    
    tfs_info <- lapply(comm_tfs$ensembl_id, function(transf){
      e_tfs <- e_in_comm %>% filter(from == transf | to == transf)
      e_reg <- tf_interactions %>% filter(tf == transf)
      
      return(list(tf = transf, non_regulatory = (nrow(e_tfs)-nrow(e_reg)), regulatory = nrow(e_reg)))
    })
    
    tfs_info <- bind_rows(tfs_info)
    tfs_info$com_id <- comm
    return(tfs_info)
  })
  all_tfs <- bind_rows(all_tfs)
  all_tfs <- all_tfs %>% left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                         by = c("tf" = "ensemblID"))
  all_tfs <- pivot_longer(all_tfs, c("non_regulatory", "regulatory"), names_to = "interaction")  
  
  comm_info <- read_tsv(paste0("data/communities/luma-communities-info-", algrthm, ".tsv")) %>%
    left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                  by = c("pg_gene" = "ensemblID"))
  all_tfs <- all_tfs %>% left_join(comm_info %>% 
                          dplyr::select(com_id, symbol), by = "com_id", 
                        suffix = c("_gene", "_community"))
  
  colors <- c("white", "#336600")
  
  p <- ggplot(all_tfs, aes(fill=interaction, y=value, x=symbol_gene)) + 
    geom_bar(position="stack", stat="identity", width = 0.85,  color = "black") +
    ylab("Interactions in network") +
    xlab("") +
    theme_few(base_size = 8) +
    theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 8), 
          legend.position="bottom", legend.key.size = unit(0.6, "line"),
          legend.title = element_blank(),
          legend.direction = "vertical") + 
    facet_grid( ~ symbol_community, scales = "free_x",  space = "free", shrink = FALSE) +
    scale_fill_manual(name = "Interaction type", breaks = c("regulatory"), 
                      values = colors, labels = c("Regulatory interactions")) 
  
  png(paste0("figures/tfs/regulatory_interactions-", algrthm, ".png"), 
      units="in", width=5, height=2.5, res=300)
  print(p)
  dev.off()
})
