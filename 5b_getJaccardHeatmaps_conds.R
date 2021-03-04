library(readr)
library(tidyr)
library(dplyr)
library(janitor)
library(ComplexHeatmap)
library(circlize)

conds <- c("healthy", "luma")
algorithms <- c("fast_greedy", "infomap" ,"leading_eigenvector", "multi_level")
labels_alg <- c("Fast Greedy", "Infomap", "Leading Eigenvector", "Louvain")
names(labels_alg) <- algorithms

biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt") %>% clean_names() %>%
  select(gene_stable_id, hgnc_symbol)

col_fun = colorRamp2(c(0, 1), c("white", "steelblue4"))

all_enrichments <- lapply(algorithms, function(alg){
  
    jaccs <- read_tsv(paste0("data/enrich/luma-healthy-", alg, ".tsv"))
    healthy_info <- read_tsv(paste0("data/communities/healthy-communities-info-", 
                                  alg, ".tsv"))
    luma_info <- read_tsv(paste0("data/communities/luma-communities-info-", 
                                  alg, ".tsv"))
    
    jaccard_matrix <- jaccs %>% 
      left_join(healthy_info %>% select(com_id, pg_gene, order),
                by=c("comm1" ="com_id")) %>%
      left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
      select(-pg_gene) %>% rename(gene_healthy = hgnc_symbol, order_healthy = order) %>% 
      left_join(luma_info %>% select(com_id, pg_gene, order),
                by=c("comm2" ="com_id")) %>%
      left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
      select(-pg_gene) %>% rename(gene_luma = hgnc_symbol, order_luma = order ) %>% 
      mutate(c1 = paste0(gene_healthy, " (", order_healthy, "-", terms1, ")"),
             c2 = paste0(gene_luma, " (", order_luma, "-", terms2, ")")) %>%
      pivot_wider(id_cols = c1, names_from = c2, values_from = jaccard) %>%
      select(sort(colnames(.)))
    
    comm1_names <- jaccard_matrix %>% select(c1) %>% unlist(use.names = F)
    jaccard_matrix <- jaccard_matrix %>%select(-c1) %>% as.matrix()
    rownames(jaccard_matrix) <- comm1_names
    jaccard_matrix <- jaccard_matrix[sort(rownames(jaccard_matrix)), ]
    
    ht <- Heatmap(jaccard_matrix, name = "Jaccard\n index", col = col_fun, 
                  cluster_rows = F, cluster_columns = F, width = ncol(jaccard_matrix), 
                  height = nrow(jaccard_matrix), 
                  column_title = paste0(labels_alg[alg], "\nLuma(", ncol(jaccard_matrix), ")"),
                  row_title = paste0("Healthy (", nrow(jaccard_matrix), ")"),
                  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                  row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12),
                  heatmap_legend_param=list(title_gp=gpar(fontsize=10), 
                                            labels_gp=gpar(fontsize=8)))
    
    png(paste0("figures/enrich/lh-", alg, "-heatmap.png"), 
        units="in", width=5.5, height=5.5, res=300)
    draw(ht)
    dev.off()  
})

