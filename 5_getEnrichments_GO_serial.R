library(readr)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

conds <- c("luma", "healthy")
algorithms <- c("leading_eigenvector", "fast_greedy", "multi_level", "infomap")

m <- mapply(function(cond, algrthm) {
  
  cat("Working with condition: ", cond, ", ", algrthm, "\n")
  
  membership <- read_tsv(file = paste0("data/communities/", cond,  "-communities-",
                                       algrthm, ".tsv"))
  gene_universe <- read_tsv(file = paste0("data/genes_in_exp_matrix.txt"), 
                            col_names = c("ensemblID"), skip = 1)
  
  gene_universe <- unlist(gene_universe$ensemblID)
  
  all_enrichments <- lapply(X = unique(membership$community),
                            FUN = function(com){
				
			cat("Working with community: ", com, "\n")

      gene_list <- unlist(membership %>% filter(community == com) %>%
                           dplyr::select(ensemblID))
      
      if(length(gene_list) >= 5) {
        
        ego <- enrichGO(gene          = gene_list,
                        universe      = gene_universe,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = "ENSEMBL",
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.005,
                        qvalueCutoff  = 0.01,
                        minGSSize     = 10,
                        readable      = FALSE)
        
        if(nrow(as.data.frame(ego)) > 0) {
         # ego <- simplify(ego)
          ego_df <- as.data.frame(ego)
          ego_df$commun <- com 
          
          ## Plots with gene Symbols 
          ## with setReadable from clusterProfile
          ## but it's not parallelizable
        
          # png(paste0("figures/enrich-universe-simple/", cond, "-go-", com, "-dotplot.png"), 
          #     width = 600, height = 800)
          # print(dotplot(ego, showCategory=20))
          # dev.off()
          #   
          # png(paste0("figures/enrich-universe-simple/", cond, "-go-", com, "-cnetplot.png"), 
          #       width = 600, height = 600)
          # print(cnetplot(ego, node_label="category", showCategory = 10))
          # dev.off()
          # 
          return(ego_df)  
        }
      }
        
      return(NULL)
  })
  
  all_enrichments <- plyr::compact(all_enrichments) 
  all_enrichments <- plyr::ldply(all_enrichments) 
  
  write.table(all_enrichments, file =  paste0("data/enrich/", cond, "-go-enrichments-",
                                              algrthm, ".tsv"),
              quote = F, row.names = F, sep = "\t")
  
}, conds, algorithms)
