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

alg_pairs <- as.data.frame(combn(algorithms, 2), stringsAsFactors = F) %>%
  bind_cols(t(data.frame(algorithms, algorithms, stringsAsFactors = F)))

all_enrichments <- lapply(conds, function(cond){
  m <- apply(alg_pairs, MARGIN = 2, FUN=function(ap)  {
    alg1 <- ap[1]
    alg2 <- ap[2]
    jaccs <- read_tsv(paste0("data/enrich/", cond, "-", alg1, "-", alg2, ".tsv"))
    comm_info1 <- read_tsv(paste0("data/communities/", cond, "-communities-info-", 
                                 alg1, ".tsv"))
    comm_info2 <- read_tsv(paste0("data/communities/", cond, "-communities-info-", 
                                  alg2, ".tsv"))
    
    jaccard_matrix <- jaccs %>% 
      left_join(comm_info1 %>% select(com_id, pg_gene),
                                          by=c("comm1" ="com_id")) %>%
      left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
      select(-pg_gene) %>% rename(gene_comm1 = hgnc_symbol) %>% 
      left_join(comm_info2 %>% select(com_id, pg_gene),
                                          by=c("comm2" ="com_id")) %>%
      left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
      select(-pg_gene) %>% rename(gene_comm2 = hgnc_symbol) %>% 
      mutate(c1 = paste0(gene_comm1, " (", terms1, ")"),
             c2 = paste0(gene_comm2, " (", terms2, ")")) %>%
      pivot_wider(id_cols = c1, names_from = c2, values_from = jaccard) %>%
      select(sort(colnames(.)))
    comm1_names <- jaccard_matrix %>% select(c1) %>% unlist(use.names = F)
    jaccard_matrix <- jaccard_matrix %>%select(-c1) %>% as.matrix()
    rownames(jaccard_matrix) <- comm1_names
    jaccard_matrix <- jaccard_matrix[sort(rownames(jaccard_matrix)), ]
    
    ht <- Heatmap(jaccard_matrix, name = "Jaccard\n index", col = col_fun, 
            cluster_rows = F, cluster_columns = F, width = ncol(jaccard_matrix), 
            height = nrow(jaccard_matrix), 
            column_title = paste0(labels_alg[alg1], " vs ", labels_alg[alg2]),
            row_names_gp = gpar(fontsize = 18), column_names_gp = gpar(fontsize = 18),
            column_title_gp = gpar(fontsize = 24))
    
    png(paste0("figures/enrich/", alg1, "-", alg2, "-heatmap.png"), 
        width = 800, height = 800)
    draw(ht, legend_title_gp = gpar(fontsize = 18), 
         legend_labels_gp = gpar(fontsize = 12))
    dev.off()  
  })
})
