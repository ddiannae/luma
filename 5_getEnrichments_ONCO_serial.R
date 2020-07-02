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

m <- lapply(conds, function(cond) {
  cat("Working with condition: ", cond, "\n")
  
  membership <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  membership$entrez <- mapIds(org.Hs.eg.db,
                              keys = membership$ensemblID,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
  
  gene_universe <- read_tsv(file = paste0("data/genes_in_exp_matrix.txt"), 
                            col_names = c("ensemblID"), skip = 1)
  gene_universe$entrez <- mapIds(org.Hs.eg.db,
                                 keys = gene_universe$ensemblID,
                                 column="ENTREZID",
                                 keytype="ENSEMBL",
                                 multiVals="first")
  gene_universe <- unlist(gene_universe$entrez)
  
  all_enrichments <- lapply(X = unique(membership$community),
                            FUN = function(com){
                                cat("Working with community: ", com, "\n")

                                gene_list <- unlist(membership %>% filter(community == com) %>%
                                                     dplyr::select(entrez))
                                
                                if(length(gene_list) >= 5) {
                                  
                                  eonco <- enricher(gene        = gene_list,
                                                  universe      = gene_universe,
                                                  pAdjustMethod = "BH",
                                                  pvalueCutoff  = 0.05,
                                                  qvalueCutoff  = 0.1,
                                                  minGSSize     = 10,
                                                  TERM2GENE = ONCO)
                                  
                                
                                  if(nrow(as.data.frame(eonco)) > 0) {
                                    eonco_df <- as.data.frame(eonco)
                                    eonco_df$commun <- com 
                                    
                                    ## Plots with gene Symbols 
                                    ## with setReadable from clusterProfile
                                    ## but it's not parallelizable
                                  
                                    png(paste0("figures/enrich-universe/", cond, "-onco-", com, "-dotplot.png"), 
                                        width = 600, height = 800)
                                    print(dotplot(eonco, showCategory=20))
                                    dev.off()
                                      
                                    png(paste0("figures/enrich-universe/", cond, "-onco-", com, "-cnetplot.png"), 
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
  
  write.table(all_enrichments, file =  paste0("data/enrich-universe/",cond, "-onco-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
  
})