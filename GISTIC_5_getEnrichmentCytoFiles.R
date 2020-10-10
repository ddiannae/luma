library(readr)
library(stringr)
library(dplyr)
library(ggplot2)

all_enrichments <- read_tsv("data/luma-gistic/community-enrichment.tsv")
gistic_peaks <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv")
gistic_peaks <- gistic_peaks %>% mutate(id = paste0(id, "_", type))
comm_info <- read_tsv("data/communities/luma-communities-info.tsv", 
                      col_types = cols_only(com_id = col_double(), order = col_double(), 
                                            pg_gene = col_character())) %>%
  left_join(read_tsv("data/network-tables/luma-20127-vertices.tsv",
                     col_types = cols_only(ensemblID = col_character(), symbol = col_character())),
            by = c("pg_gene" = "ensemblID"))

comm_expr <- read_tsv("data/assortativity/luma-exp-assortativity.tsv") %>% 
  select(community_id, totalfrac, mean_diff_exp)
comm_chr <- read_tsv("data/assortativity/luma-chr-assortativity.tsv",
                     col_types = cols_only(community_id = col_double(),
                                           totalfrac = col_double()))

luma_assort <- comm_chr %>% 
  inner_join(comm_expr, by = "community_id", suffix = c("_chr", "_exp"))

comm_info <- comm_info %>% inner_join(luma_assort, by = c( "com_id" = "community_id"))

all_lesions <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv")
all_lesions <- all_lesions %>% mutate(id = paste0(id, "_", type))

comm_info <- comm_info %>% semi_join(all_enrichments, by = c("com_id" = "comm"))
gistic_peaks <- gistic_peaks %>% semi_join(all_enrichments, by = c("id" = "ID"))

comm_info <- comm_info %>% rename(id = com_id)
comm_info$node_type <- "comm"

gistic_peaks$node_type <- "peak"
write_tsv(comm_info, "data/luma-gistic/enriched-comms.tsv")
write_tsv(gistic_peaks, "data/luma-gistic/enriched-peaks.tsv")

luma_genes <- read_tsv("data/network-tables/luma-20127-vertices.tsv")
lesion_genes <- read_tsv("data/luma-gistic/luma-network-genes.tsv")

### Get genes in enriched communities
luma_genes <- luma_genes %>% filter(community %in% comm_info$id)
luma_genes %>% select( ensemblID) %>% write_tsv("data/luma-gistic/genes-in-enriched-comms.tsv")
### Select only genes in lesions
lesion_genes <- luma_genes %>% filter(ensemblID %in% lesion_genes$ensembl_id)
lesion_genes %>% select( ensemblID) %>% write_tsv("data/luma-gistic/genes-in-enriched-comms-lesion.tsv")
