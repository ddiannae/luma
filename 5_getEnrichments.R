library(readr)
library(dplyr)
library(parallel)
library(DBI)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

conds <- c("healthy", "luma")

enrichCommunities <- function(membership, universe, filename, 
                              T2GDB, T2NDB, G2S, simplifyGO = FALSE){
  
  all_enrichments <- mclapply(X = unique(membership$community), 
                              mc.cores = 75,
                              FUN = function(com){
                                
      geneList <- unlist(membership %>% filter(community == com) %>%
                           dplyr::select(ensemblID))
      
      if(length(geneList) >= 5) {
        
        ego <- enricher(gene          = geneList,
                        universe      = universe,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        minGSSize     = 10,
                        TERM2GENE     = T2GDB,
                        TERM2NAME = T2NDB)
        
        if(nrow(as.data.frame(ego)) > 0) {
          
          if(simplifyGO == TRUE) {
            ego@ontology = "BP"
            ego@organism = "Homo sapiens"
            
            ego <- simplify(
              ego,
              cutoff = 0.7,
              by = "p.adjust",
              select_fun = min,
              measure = "Wang",
              semData = NULL
            )  
          }
          
          ego_df <- as.data.frame(ego)
          ego_df$commun <- com
          
          ## Plots with gene Symbols 
          ## with setReadable from clusterProfile
          ## but it's not parallelizable
          
          png(paste0("figures/enrich/", filename, "-", com, "-dotplot.png"), 
                     width = 600, height = 800)
          print(dotplot(ego, showCategory=20))
          dev.off()
          
          gs <- G2S[G2S[[1]] %in% ego@gene, ]
          
          ego@gene2Symbol <- setNames(gs[[2]], gs[[1]])
          ego@readable = TRUE
          ego@keytype = "ENSEMBL"
          
          res <- ego@result
          gn <- ego@gene2Symbol
          gc <- geneInCategory(ego)
          gc <- lapply(gc, function(i) gn[i])
          gc <- gc[as.character(res$ID)]
          geneID <- sapply(gc, paste0, collapse="/")
          res$geneID <- unlist(geneID)
          
          ego@result <- res
         
          png(paste0("figures/enrich/", filename, "-", com, "-cnetplot.png"), 
              width = 600, height = 600)
          print(cnetplot(ego, node_label="all"))
          dev.off()
          
          return(ego_df)  
        }
      }
      return(NULL)
    })
  
  all_enrichments <- plyr::compact(all_enrichments) 
  all_enrichments <- plyr::ldply(all_enrichments) 
  
  write.table(all_enrichments, file =  paste0("data/enrich/", filename,  "-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
}

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

##### GET GO DB #####
GO_BP <- getGO_DB()
GO_BP <- tidyr::unnest(tibble::enframe(GO_BP), cols = c(value))
GO_BP_names <- go2term(GO_BP$name)

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
  gene_sets <- read_tsv(file = paste0("data/", cond,  "-communities.tsv"))
  gene_sets$symbol <- mapIds(org.Hs.eg.db,
                             keys = gene_sets$ensemblID,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
  
  comm_info <- read_tsv(file = paste0("data/", cond,  "-communities-info.tsv"))
  comm_info <- comm_info[, c("com_id", "pg_gene")]
  comm_info$symbol <- mapIds(org.Hs.eg.db,
                            keys = comm_info$pg_gene,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
  
  gene_universe <- read_tsv(file = paste0("data/network-tables/", cond,  "-vertices.tsv"))
  gene_universe <- unlist(gene_universe$ensemblID)

  enrichCommunities(membership = gene_sets[, c("ensemblID", "community")],
                    universe = gene_universe,
                    filename = paste0(cond, "-kegg"),
                    T2GDB = KEGG_db[, c("from", "ensembl")], 
                    T2NDB = KEGG$KEGGPATHID2NAME,
                    G2S = gene_sets[, c("ensemblID", "symbol")])
  
  enrichCommunities(membership = gene_sets[, c("ensemblID", "community")],
                    universe = gene_universe,
                    filename = paste0(cond, "-go"),
                    T2GDB = GO_BP, 
                    T2NDB = GO_BP_names,
                    G2S = gene_sets[, c("ensemblID", "symbol")])
  
  enrichCommunities(membership = gene_sets[, c("ensemblID", "community")],
                    universe = gene_universe,
                    filename = paste0(cond, "-onco"),
                    T2GDB = ONCO[, c("gs_name", "ensembl")], 
                    T2NDB = NULL,
                    G2S = gene_sets[, c("ensemblID", "symbol")])
  
})
