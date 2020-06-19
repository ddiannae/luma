library(readr)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

conds <- c("healthy", "luma")

#### GET KEGG DB ###
KEGG <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
KEGG_db <- KEGG$KEGGPATHID2EXTID

KEGG_db$ensembl <- mapIds(org.Hs.eg.db,
                          keys = KEGG_db$to,
                          column="ENSEMBL",
                          keytype="ENTREZID",
                          multiVals="list")

KEGG_db <- tidyr::unnest(as_tibble(KEGG_db), cols = c(ensembl))

### GET ONCO SETS DB ###
ONCO <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
ONCO$ensembl <- mapIds(org.Hs.eg.db,
                       keys = as.character(ONCO$entrez_gene),
                       column="ENSEMBL",
                       keytype="ENTREZID",
                       multiVals="list")
ONCO <- tidyr::unnest(ONCO, cols = c(ensembl))

m <- lapply(conds, function(cond) {
  membership <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  gene_universe <- read_tsv(file = paste0("data/network-tables/", cond,  "-vertices.tsv"))
  gene_universe <- unlist(gene_universe$ensemblID)
  
  all_enrichments <- lapply(X = unique(membership$community),
                            FUN = function(com){

                                gene_list <- unlist(membership %>% filter(community == com) %>%
                                                     dplyr::select(ensemblID))
                                
                                if(length(gene_list) >= 5) {
                                  
                                  ego <- enrichGO(gene          = gene_list,
                                                  universe      = gene_universe,
                                                  OrgDb         = org.Hs.eg.db,
                                                  keyType       = "ENSEMBL",
                                                  ont           = "BP",
                                                  pAdjustMethod = "BH",
                                                  pvalueCutoff  = 0.01,
                                                  qvalueCutoff  = 0.05,
                                                  minGSSize     = 10,
                                                  readable      = FALSE)
                                  
                                
                                  if(nrow(as.data.frame(ego)) > 0) {
                                    ego_df <- as.data.frame(ego)
                                    ego_df$commun <- com 
                                    
                                    ## Plots with gene Symbols 
                                    ## with setReadable from clusterProfile
                                    ## but it's not parallelizable
                                  
                                    png(paste0("figures/enrich-serial/", filename, "-", com, "-dotplot.png"), 
                                        width = 600, height = 800)
                                    print(dotplot(ego, showCategory=20))
                                    dev.off()
                                      
                                    png(paste0("figures/enrich-serial/", filename, "-", com, "-cnetplot.png"), 
                                          width = 600, height = 600)
                                    print(cnetplot(ego, node_label="category", showCategory = 10))
                                    dev.off()
                                    
                                    return(ego_df)  
                                  }
                                }
                                  
                                return(NULL)
                              })
  
  all_enrichments <- plyr::compact(all_enrichments) 
  all_enrichments <- plyr::ldply(all_enrichments) 
  
  write.table(all_enrichments, file =  paste0("data/enrich-serial/",cond, "go-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
  
  
})