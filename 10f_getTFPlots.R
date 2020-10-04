library(igraph)
library(readr)
library(ggplot2)
library(tidyr)
library(ggthemes)

communities <- read_tsv("data/enrich-universe/communities-more-30terms.tsv") %>% unlist()
interactions <- read_tsv("data/network-tables/luma-20127-interactions.tsv")
vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
colnames(interactions)[1] = "from"
colnames(interactions)[2] = "to"
tf_interactions <- read_tsv("data/tfs/gtrd/all_tfs_interactions.tsv")
biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt",   
          col_names = c("ensemblID", "chr", "start", "end",  "gc", "type", "symbol"), skip = 1)

all_tfs <- lapply(communities, function(comm){
  v_in_comm <- vertices %>% filter(community == comm) %>% 
    dplyr::select(ensemblID) %>% unlist(use.names = F)
  comm_tfs <- read_tsv(paste0("data/tfs/gtrd/tfs_in_comm_", comm, ".txt"))
  e_in_comm <- interactions %>% filter(from %in% v_in_comm & to %in% v_in_comm)
  
  tfs_info <- lapply(comm_tfs$ensembl_id, function(tf){
    e_tfs <- e_in_comm %>% filter(from == tf | to == tf)
    e_reg <- tf_interactions %>% filter(from == tf | to == tf)
    
    return(list(tf = tf, non_regulatory = (nrow(e_tfs)-nrow(e_reg)), regulatory = nrow(e_reg)))
  })
  
  tfs_info <- bind_rows(tfs_info)
  tfs_info$com_id <- comm
  return(tfs_info)
})
all_tfs <- bind_rows(all_tfs)
all_tfs <- all_tfs %>% left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                       by = c("tf" = "ensemblID"))
all_tfs <- pivot_longer(all_tfs, c("non_regulatory", "regulatory"), names_to = "interaction")  

comm_info <- read_tsv("data/communities/luma-communities-info.tsv") %>%
  left_join(biomart %>% dplyr::select(ensemblID, symbol), 
                by = c("pg_gene" = "ensemblID"))
all_tfs <- all_tfs %>% left_join(comm_info %>% 
                        dplyr::select(com_id, symbol), by = "com_id", 
                      suffix = c("_gene", "_community"))

colors <- c("white", "darkolivegreen4")

png(paste0("figures/tfs/regulatory_interactions.png"), width= 1150, height = 600)
ggplot(all_tfs, aes(fill=interaction, y=value, x=symbol_gene)) + 
  geom_bar(position="stack", stat="identity", width = .8,  color = "black") +
  ylab("Interactions in network") +
  xlab("") +
  theme_few(base_size = 20) +
  theme(axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.direction = "vertical") + 
  facet_grid( ~ symbol_community, scales = "free_x",  space = "free", shrink = FALSE) +
  scale_fill_manual(name = "Interaction type", breaks = c("regulatory"), 
                    values = colors, labels = c("Regulatory interactions")) 
dev.off()
