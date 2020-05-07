library(readr)
library(dplyr)
library(parallel)
library(DBI)
library(org.Hs.eg.db)
library(clusterProfiler)

conds <- c("healthy", "luma")

## https://github.com/YuLab-SMU/clusterProfiler/issues/207
## Creamos el objeto GO_BP porque: https://support.bioconductor.org/p/38541/

## esta funcion la adaptamos de la funcion annFUN.org del paquete topGO
## no podemos usar topGO porque no puede paralelizarse
getGO_DB <- function() {
  .sql <- paste("SELECT DISTINCT ensembl_id, go_id FROM ensembl",
                " INNER JOIN go_bp", 
                " USING(_id)", sep = "")
  retVal <- dbGetQuery(get("org.Hs.eg_dbconn")(), .sql)
  
  ## split the table into a named list of GOs
  return(split(retVal[["ensembl_id"]], retVal[["go_id"]]))
}

GO_BP <- getGO_DB()
GO_BP <- tibble::enframe(GO_BP)

m <- lapply(conds, function(cond) {
  gene_sets <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  gene_universe <- read_tsv(file = paste0("data/network-tables/", cond,  "-vertices.tsv"))
  gene_universe <- unlist(gene_universe$ensemblID)

  all_enrichments <- mclapply(X = unique(gene_sets$community), 
                              mc.cores = 75,
                              FUN = function(com){
   
      geneList <- unlist(gene_sets %>% filter( community == com) %>%
                           dplyr::select(ensemblID))
      
      ego <- enricher(gene          = geneList,
                      universe      = gene_universe,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      minGSSize     = 10,
                      TERM2GENE     = GO_BP)
      
      if(nrow(as.data.frame(ego)) > 0) {
        ego@ontology = "BP"
        ego@organism = "Homo sapiens"
        
        ego <- as.data.frame(simplify(
          ego,
          cutoff = 0.7,
          by = "p.adjust",
          select_fun = min,
          measure = "Wang",
          semData = NULL
        ))
        
        ego <- ego %>% inner_join(go2term(ego$ID), by = c("ID" = "go_id")) %>% 
                dplyr::select(-Description) %>% dplyr::select(ID, Term, everything())
        
        
        ego$commun <- com
        return(ego)  
      }
      return(NULL)
  })
  
  all_enrichments <- plyr::compact(all_enrichments) 
  all_enrichments <- plyr::ldply(all_enrichments) 
 
  write.table(all_enrichments, file =  paste0("data/", cond,  "-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
})

