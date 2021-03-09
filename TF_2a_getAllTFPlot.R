library(igraph)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)


biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt",   
                    col_names = c("ensemblID", "chr", "start", "end",  "gc", "type", "symbol"), skip = 1)

algrthm <- "multi_level"
conds <- c("healthy", "luma")

m <- lapply(conds, function(cond) {
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-20127-interactions.tsv"))
  vertices <- read_tsv(paste0("data/network-tables/", cond, "-20127-vertices.tsv"))
  
  colnames(interactions)[1] = "from"
  colnames(interactions)[2] = "to"
  
  gene_comm <- read_tsv(paste0("data/communities/", cond, "-communities-", algrthm, ".tsv"))
  communities <- unique(gene_comm$community)
  tf_interactions <- read_tsv(paste0("data/tfs/gtrd/", cond, "-tfs_interactions-", algrthm, ".tsv"))
  
  all_tfs <- lapply(communities, function(comm){
    
    v_in_comm <- gene_comm %>% filter(community == comm) %>% 
      dplyr::select(ensemblID) %>% unlist(use.names = F)
    
    comm_tfs <- read_tsv(paste0("data/tfs/gtrd/", cond, "-tfs_in_comm_", comm, "-", algrthm, ".txt"))
    e_in_comm <- interactions %>% filter(from %in% v_in_comm & to %in% v_in_comm)

      
      tfs_info <- lapply(comm_tfs$ensembl_id, function(transf){
        e_tfs <- e_in_comm %>% filter(from == transf | to == transf)
        e_reg <- tf_interactions %>% filter(tf == transf)
        
        return(list(tf = transf, non_regulatory = (nrow(e_tfs)-nrow(e_reg)), regulatory = nrow(e_reg)))
      })
      
      tfs_info <- bind_rows(tfs_info)
    if(nrow(tfs_info) > 0 ){
      tfs_info$com_id <- comm
      tfs_info <- tfs_info %>% left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                                         by = c("tf" = "ensemblID"))
      tfs_info <- pivot_longer(tfs_info, c("non_regulatory", "regulatory"), names_to = "interaction")  
      
      comm_info <- read_tsv(paste0("data/communities/luma-communities-info-", algrthm, ".tsv")) %>%
        left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                  by = c("pg_gene" = "ensemblID"))
      tfs_info <- tfs_info %>% left_join(comm_info %>% 
                                           dplyr::select(com_id, symbol), by = "com_id", 
                                         suffix = c("_gene", "_community"))
      colors <- c("white", "#336600")
      
      p <- ggplot(tfs_info, aes(fill=interaction, y=value, x=symbol_gene)) + 
        geom_bar(position="stack", stat="identity", width = 0.85,  color = "black") +
        ylab("Interactions in network") +
        xlab("") +
        theme_few(base_size = 4) +
        theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position="bottom", legend.key.size = unit(0.6, "line"),
              legend.title = element_blank(),
              legend.direction = "vertical") + 
        scale_fill_manual(name = "Interaction type", breaks = c("regulatory"), 
                          values = colors, labels = c("Regulatory interactions")) 
      
      png(paste0("figures/tfs/", cond, "-regulatory_interactions-", comm, "-", algrthm, ".png"), 
          units="in", width=3, height=1.5, res=300)
      print(p)
      dev.off()
    }
    
  })
})
