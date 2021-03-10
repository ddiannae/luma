library(readr)
library(tidyr)
library(dplyr)
library(tidyr)
library(igraph)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(janitor)

algrthm = "multi_level"
conds <- c("healthy", "luma")
biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt") %>% clean_names() %>%
  select(gene_stable_id, hgnc_symbol) %>% rename("ensemblID" = "gene_stable_id", "symbol" = "hgnc_symbol")


all_inter <- lapply(conds, function(cond){
  gene_comm <- read_tsv(paste0("data/tfs/gtrd/", cond , "-tfs_interactions-", algrthm, ".tsv"))
  
  tf_inter <- lapply(unique(gene_comm$community), function(comm){
    comm_tf_inter <- read_tsv(paste0("data/tfs/gtrd/", cond, "-all_tfs_interactions-", comm, ".tsv"))
  })
  tf_inter <- bind_rows(tf_inter)
  return(tf_inter)
})
all_inter <- bind_rows(all_inter)
all_inter <- distinct(all_inter )
all_inter  <- all_inter  %>% left_join(biomart, by = c("from" ="ensemblID")) %>% 
  select(symbol, to, from)


all_enr <-  lapply(conds, function(cond) {
  
  interactions <- read_tsv(paste0("data/network-tables/", cond, "-20127-interactions.tsv"))
  vertices <- read_tsv(paste0("data/network-tables/", cond, "-20127-vertices.tsv"))
  colnames(interactions)[1:2] <- c("from", "to")
  universe <- vertices$ensemblID
  net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)
  tfs_in_net <-  vertices %>% filter(ensemblID %in% unique(all_inter$from)) %>%
    select(ensemblID) %>% unlist()
  enrichs <-  parallel::mclapply(tfs_in_net,  mc.cores = 70, FUN = function(tf) {
      
      geneList <- unlist(neighbors(net, tf)$name)
      etf <- enricher(gene  = geneList,
                          universe  = universe,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.1,
                          minGSSize     = 1,
                          maxGSSize = 10000,
                          TERM2GENE     = all_inter)
      etf <- as.data.frame(etf)
      
      if(nrow(etf) > 0) {
        etf$tf <- tf
        return(etf) 
      }
      return(NULL)
  })
    enrichs <- bind_rows(enrichs)
   
    if(nrow(enrichs) > 0) {
      enrichs$cond <- cond
      return(enrichs) 
    }
    return(NULL)
})
all_enr <- bind_rows(all_enr)
write_tsv(all_enr, file = "data/tfs/gtrd/all-tf-enrichments.tsv")


### ------------------------------------------------
all_enr <- read_tsv("data/tfs/gtrd/all-tf-enrichments.tsv")

healthy_enr <- all_enr %>% filter(cond == "healthy")
luma_enr <- all_enr %>% filter(cond == "luma")

luma_enr <- luma_enr %>%
  left_join(read_tsv("data/network-tables/luma-20127-vertices.tsv",
                     col_types = cols_only(ensemblID = col_character(), symbol= col_character())),
            by = c("tf" = "ensemblID"))
luma_enr_matrix <- luma_enr %>%
  pivot_wider(id_cols = ID, names_from = symbol, values_from = p.adjust, values_fill = 0)

tfs_bd <- luma_enr_matrix[,1] %>% unlist()
luma_enr_matrix <- luma_enr_matrix %>% select(-ID) %>% as.matrix()
rownames(luma_enr_matrix) <- tfs_bd

col_fun = colorRamp2(c(0, 0 +1e-16, 0.05) , c("white", "steelblue4",  "white"))
ht <- Heatmap(luma_enr_matrix, name = "p-value", col = col_fun,
              row_title = paste0("TFs in data base"),
              column_title = paste0("TFs in LumA GCN"),
              row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 7)
              )
png(paste0("figures/tfs/luma-tfs-heatmap.png"),
    units="in", width=5.5, height=12, res=300)
draw(ht)
dev.off()

healthy_enr <- healthy_enr %>%
  left_join(read_tsv("data/network-tables/healthy-20127-vertices.tsv",
                     col_types = cols_only(ensemblID = col_character(), symbol= col_character())),
            by = c("tf" = "ensemblID"))
healthy_enr_matrix <- healthy_enr %>%
  pivot_wider(id_cols = ID, names_from = symbol, values_from = p.adjust, values_fill = 0)

tfs_bd <- healthy_enr_matrix[,1] %>% unlist()
healthy_enr_matrix <- healthy_enr_matrix %>% select(-ID) %>% as.matrix()
rownames(healthy_enr_matrix) <- tfs_bd

col_fun = colorRamp2(c(0, 0 +1e-16, 0.05) , c("white", "steelblue4",  "white"))
ht <- Heatmap(healthy_enr_matrix,  name = "p-value", col = col_fun,
              row_title = paste0("TFs in data base"),
              column_title = paste0("TFs in Healthy GCN"),
              row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12),
              row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7)
)
png(paste0("figures/tfs/healthy-tfs-heatmap.png"),
    units="in", width=5, height=8, res=300)
draw(ht)
dev.off()

