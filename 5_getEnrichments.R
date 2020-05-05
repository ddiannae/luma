library(topGO)
library(dplyr)

conds <- c("healthy", "luma")

gene_list <- lapply(conds, function(cond) {
  genes <- read_tsv(file = paste0("data/network-tables/", cond,  "-vertices.tsv"))
  return(genes$ensemblID)
})
gene_list <- unique(unlist(gene_list))

m <- lapply(conds, function(cond) {
  gene_sets <- read_tsv(file = paste0("data/network-tables/", cond,  "-communities.tsv"))
  
  all_enrichments <- lapply(unique(gene_sets$community), function(com){
   
      geneList <- unlist(gene_sets %>% filter( community == com) %>%
                           dplyr::select(ensemblID))
      geneList <- setNames(factor(as.integer(gene_list %in% geneList)), gene_list)
      
      GOdata <- new("topGOdata",
                    description = "GO analysis of genes in community",
                    ontology = "BP",
                    allGenes = geneList,
                    nodeSize = 10,
                    annot = annFUN.org,
                    ID = "Ensembl", 
                    mapping = "org.Hs.eg")
      
      resultFisher <- runTest(GOdata, "elim","fisher")
      ss <- score(resultFisher)
      cat("Community ", com, " has ", sum(ss < 0.001), " significant processes")
      resultsTable <- GenTable(GOdata, classicFisher = resultFisher, numChar = 200,
                               topNodes = sum(ss < 0.001))
      resultsTable$commun <- com
      return(resultsTable)
  })
  all_enrichments <- bind_rows(all_enrichments)
 
  write.table(all_enrichments, file =  paste0("data/", cond,  "-enrichments.tsv"),
              quote = F, row.names = F, sep = "\t")
})


