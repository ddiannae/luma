library(msigdbr)
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DBI)

luma_deg <- read_tsv("data/luma-deg-ebayes.tsv")
luma_deg <- luma_deg %>% filter(FDR < 0.05)

geneList <- luma_deg$coef
names(geneList) <- luma_deg$gene
geneList <- sort(geneList, decreasing = TRUE)

###### ONCOGENIC GENE SETS
##  https://www.gsea-msigdb.org/gsea/msigdb
##
######

ONCO <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
ONCO$ensembl <- mapIds(org.Hs.eg.db,
                          keys = as.character(ONCO$entrez_gene),
                          column="ENSEMBL",
                          keytype="ENTREZID",
                          multiVals="list")
ONCO <- tidyr::unnest(ONCO, cols = c(ensembl))

gsea_onco <- GSEA(geneList, 
                  TERM2GENE = ONCO[, c("gs_name", "ensembl")],
                  nPerm = 10000, 
                  pvalueCutoff = 0.05)

png("figures/onco_dotplot.png", width = 800, height = 1000)
print(dotplot(gsea_onco, showCategory=20))
dev.off()

png("figures/onco_ridgeplot.png", width = 1200, height = 1200)
print(ridgeplot(gsea_onco))
dev.off()

png("figures/onco_gseaplot.png", width = 1500, height = 800)
print(gseaplot2(gsea_onco, geneSetID = 1:5))
dev.off()

write_tsv(as.data.frame(gsea_onco), path = "data/GSEA/oncogenic_genesets_gsea.tsv")

###### GENE ONTOLOGY BIO PROCESSES
##  
##
######

getGO_DB <- function() {
  .sql <- paste("SELECT DISTINCT ensembl_id, go_id FROM ensembl",
                " INNER JOIN go_bp", 
                " USING(_id)", sep = "")
  retVal <- dbGetQuery(get("org.Hs.eg_dbconn")(), .sql)
  
  ## split the table into a named list of GOs
  return(split(retVal[["ensembl_id"]], retVal[["go_id"]]))
}

GO_BP <- getGO_DB()
GO_BP <- tidyr::unnest(tibble::enframe(GO_BP), cols = c(value))
GO_BP_names <- go2term(GO_BP$name)

gsea_go <- GSEA(geneList, TERM2GENE = GO_BP, 
                TERM2NAME = GO_BP_names,
                nPerm = 10000, 
                pvalueCutoff = 0.05)

png("figures/go_dotplot.png", width = 800, height = 1000)
print(dotplot(gsea_go, showCategory=20))
dev.off()

png("figures/go_ridgeplot.png", width = 1200, height = 1200)
print(ridgeplot(gsea_go))
dev.off()

png("figures/go_gseaplot.png", width = 1500, height = 800)
print(gseaplot2(gsea_go, geneSetID = 1:5))
dev.off()

write_tsv(as.data.frame(gsea_go), path = "data/GSEA/go_gsea.tsv")

###### KEGG
##  
##
######

KEGG <- download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")
KEGG_db <- KEGG$KEGGPATHID2EXTID

KEGG_db$ensembl <- mapIds(org.Hs.eg.db,
                          keys = KEGG_db$to,
                          column="ENSEMBL",
                          keytype="ENTREZID",
                          multiVals="list")

KEGG_db <- tidyr::unnest(as_tibble(KEGG_db), cols = c(ensembl))

gsea_kegg <- GSEA(geneList,
                TERM2GENE = KEGG_db[, c("from", "ensembl")],
                TERM2NAME = KEGG$KEGGPATHID2NAME,
                nPerm = 10000, 
                pvalueCutoff = 0.05)

png("figures/kegg_dotplot.png", width = 800, height = 1000)
print(dotplot(gsea_kegg, showCategory=20))
dev.off()

png("figures/kegg_ridgeplot.png", width = 1200, height = 1200)
print(ridgeplot(gsea_kegg))
dev.off()

png("figures/kegg_gseaplot.png", width = 1500, height = 800)
print(gseaplot2(gsea_kegg, geneSetID = 1:5))
dev.off()

write_tsv(as.data.frame(gsea_kegg), path = "data/GSEA/kegg_gsea.tsv")
