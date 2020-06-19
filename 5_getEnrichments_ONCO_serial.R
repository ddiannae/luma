library(readr)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

conds <- c("healthy", "luma")

### GET ONCO SETS DB ###
ONCO <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)

ONCO$ensembl <- mapIds(org.Hs.eg.db,
                       keys = as.character(ONCO$entrez_gene),
                       column="ENSEMBL",
                       keytype="ENTREZID",
                       multiVals="list")
ONCO <- tidyr::unnest(ONCO, cols = c(ensembl))
ONCO <- ONCO %>% dplyr::select(gs_name, ensembl)

m <- lapply(conds, function(cond) {
  membership <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  
  gene_universe <- read_tsv(file = paste0("data/network-tables/", cond,  "-vertices.tsv"))
  gene_universe <- unlist(gene_universe$ensemblID)
  
  all_enrichments <- lapply(X = unique(membership$community),
                            FUN = function(com){

                                gene_list <- unlist(membership %>% filter(community == com) %>%
                                                     dplyr::select(ensemblID))
                                
                                if(length(gene_list) >= 5) {
                                  
                                  eonco <- enricher(gene        = gene_list,
                                                  universe      = gene_universe,
                                                  pAdjustMethod = "BH",
                                                  pvalueCutoff  = 0.01,
                                                  qvalueCutoff  = 0.05,
                                                  minGSSize     = 10,
                                                  TERM2GENE = ONCO)
                                  
                                
                                  if(nrow(as.data.frame(eonco)) > 0) {
                                    eonco_df <- as.data.frame(eonco)
                                    eonco_df$commun <- com 
                                    
                                    ## Plots with gene Symbols 
                                    ## with setReadable from clusterProfile
                                    ## but it's not parallelizable
                                  
                                    png(paste0("figures/enrich-serial/onco-", com, "-dotplot.png"), 
                                        width = 600, height = 800)
                                    print(dotplot(eonco, showCategory=20))
                                    dev.off()
                                      
                                    png(paste0("figures/enrich-serial/onco-", com, "-cnetplot.png"), 
                                          width = 600, height = 600)
                                    print(cnetplot(eonco, node_label="category", showCategory = 10))
                                    dev.off()
                                    
                                    return(eonco_df)  
                                  }
                                }
                                  
                                return(NULL)
                              })
  
  all_enrichments <- plyr::compact(all_enrichments) 
  all_enrichments <- plyr::ldply(all_enrichments) 
  
  write.table(all_enrichments, file =  paste0("data/enrich-serial/",cond, "-onco-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
  
})