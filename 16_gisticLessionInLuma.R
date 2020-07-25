library(readr)
library(dplyr)

luma_genes <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
genes_gistic <- read_tsv("data/luma-gistic/all-genes-in-lesions.tsv") %>% 
  rename(lesion_id = id , lesion_type = type)
lesions_gistic <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv")
colnames(lesions_gistic) <- paste0("lesion_",colnames(lesions_gistic))

genes_gistic <- genes_gistic %>% left_join(lesions_gistic, by = c("lesion_id", "lesion_type"))
genes_gistic <- genes_gistic %>% semi_join(luma_genes, by = c("ensembl_id" = "ensemblID"))

genes_gistic %>% select(ensembl_id, lesion_id,  lesion_q, lesion_residual_q, lesion_type_log10_q,
                        lesion_type, lesion_wide_peak ) %>% 
  write_tsv("data/luma-gistic/luma-network-genes.tsv")
c