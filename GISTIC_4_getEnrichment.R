library(readr)
library(dplyr)
library(clusterProfiler)

algrthm = "multi_level"

luma_genes <- read_tsv("data/network-tables/luma-20127-vertices.tsv") %>% 
  select(ensemblID, symbol, community)

genes_gistic <- read_tsv("data/luma-gistic/all-genes-in-lesions.tsv") %>% 
  rename(lesion_id = id , lesion_type = type)

lesions_gistic <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv")
colnames(lesions_gistic) <- paste0("lesion_",colnames(lesions_gistic))

lesion2gene <- genes_gistic %>% mutate(lesion_id = paste0(lesion_id, "_", lesion_type)) %>%
  select(lesion_id, ensembl_id)

lesion2name <- lesions_gistic %>% mutate(lesion_id = paste0(lesion_id, "_", lesion_type), 
                                         name = paste0(lesion_cytoband, "_", lesion_type)) %>%
  select(lesion_id, name)

lesion2gene %>% group_by(lesion_id) %>% tally() %>% arrange(n)
universe <- luma_genes$ensemblID

gene_comm <- read_tsv(paste0("data/communities/luma-communities-", algrthm, ".tsv"))

comms <- unique(gene_comm$community)

all_comms <- lapply(comms, function(comm) {
  
  geneList <- unlist(gene_comm %>% filter(community == comm) %>%
                       dplyr::select(ensemblID))
  
  egistic <- enricher(gene  = geneList,
                  universe  = universe,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  minGSSize     = 1,
                  TERM2GENE     = lesion2gene,
                  TERM2NAME     = lesion2name)
  egistic <- as.data.frame(egistic)
  
  if(nrow(egistic)) {
    egistic$comm <- comm
    return(egistic)
  }
})
all_comms <- bind_rows(all_comms)
all_comms %>% select(comm) %>% unique() %>% count()
write_tsv(all_comms, "data/luma-gistic/community-enrichment.tsv") 
