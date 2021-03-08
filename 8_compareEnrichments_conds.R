library(readr)
library(dplyr)
# 
# healthy_kegg <- read_tsv("./data/enrich-universe/healthy-kegg-enrichments.tsv")
# luma_kegg <- read_tsv("./data/enrich-universe/luma-kegg-enrichments.tsv")
# 
# healthy_kegg <- healthy_kegg %>% filter(p.adjust < 0.005)
# luma_kegg <- luma_kegg %>% filter(p.adjust < 0.005)
# 
# healthy_kegg <- healthy_kegg %>% select(ID, Description) %>% unique()
# luma_kegg <- luma_kegg %>% select(ID, Description) %>% unique()
# 
# shared <- healthy_kegg %>% semi_join(luma_kegg)
# luma_only <- luma_kegg %>% anti_join(healthy_kegg)
# healthy_only <- healthy_kegg %>% anti_join(luma_kegg)
# 
# write_tsv(luma_only %>% select(ID), "data/enrich-universe/luma-only-KEGG.tsv")
# write_tsv(healthy_only %>% select(ID), "data/enrich-universe/healthy-only-KEGG.tsv")
# write_tsv(shared %>% select(ID), "data/enrich-universe/shared-only-KEGG.tsv")

algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector", "multi_level")

m <- lapply(algorithms, function(algrthm) {
  
  ### Data for alluvial
  luma_enr <- read_tsv(paste0("data/enrich/luma-go-enrichments-", algrthm, ".tsv"))
  luma_intra_comm <- read_tsv(paste0("data/communities/luma-intra-communities-", algrthm, ".tsv"))
  
  healthy_enr <- read_tsv(paste0("data/enrich/healthy-go-enrichments-", algrthm, ".tsv"))
  healthy_intra_comm <- read_tsv(paste0("data/communities/healthy-intra-communities-", algrthm, ".tsv"))
  
  
  healthy_enr_min <- healthy_enr %>% select(ID, Description) %>% unique()
  luma_enr_min <- luma_enr %>% select(ID, Description) %>% unique()
  
  shared <- healthy_enr_min %>% semi_join(luma_enr_min)
  luma_only <- luma_enr_min %>% anti_join(healthy_enr_min)
  healthy_only <- healthy_enr_min %>% anti_join(luma_enr_min)
  
  write_tsv(luma_only %>% select(ID), paste0("data/enrich/luma-only-GO-", algrthm, ".tsv"))
  write_tsv(healthy_only %>% select(ID), paste0("data/enrich/healthy-only-GO-", algrthm, ".tsv"))
  write_tsv(shared %>% select(ID), paste0("data/enrich/shared-only-GO-", algrthm, ".tsv"))
  
  luma_enr %>% filter(ID %in% luma_only$ID) %>% 
    write_tsv(paste0("data/enrich/luma-only-all-info-GO-", algrthm, ".tsv"))
  healthy_enr %>% filter(ID %in% healthy_only$ID) %>% 
    write_tsv(paste0("data/enrich/healthy-only-all-info-GO-", algrthm, ".tsv"))
  healthy_enr %>% filter(ID %in% shared$ID) %>% 
    write_tsv(paste0("data/enrich/shared-only-all-info-GO-", algrthm, ".tsv"))
  
  
  luma_enr <- luma_enr %>% mutate(community_type = 
                              ifelse(commun %in% luma_intra_comm$community_id, 
                                "*", "-"))
  intra_enrich <- luma_enr %>% filter(community_type == "*")
  ## 136 términos en 9 comunidades cis
  inter_enrich <-  luma_enr %>% filter(community_type == "-")
  ## 792 términos en 20 comunidades trans
  
  healthy_enr <- healthy_enr %>% mutate(community_type = 
                            ifelse(commun %in% healthy_intra_comm$community_id, 
                                           "*", "-"))
  
  healthy_comm_info <- read_tsv(paste0("data/communities/healthy-communities-info-", algrthm, ".tsv"))
  healthy_vertices <- read_tsv("data/network-tables/healthy-20127-vertices.tsv", 
                            col_types = cols_only(ensemblID = col_character(), symbol = col_character()))
  
  luma_comm_info <- read_tsv(paste0("data/communities/luma-communities-info-", algrthm, ".tsv"))
  luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                            col_types = cols_only(ensemblID = col_character(), symbol = col_character()))
  
  
  luma_comm_info <- luma_comm_info %>% left_join(luma_vertices, by = c("pg_gene" = "ensemblID"))
  healthy_comm_info <- healthy_comm_info %>% left_join(healthy_vertices, by = c("pg_gene" = "ensemblID"))
  
  luma_enr %>% left_join(luma_comm_info %>% select(com_id, symbol), by = c("commun" = "com_id")) %>% 
    select(ID, symbol, community_type)%>% unique() %>% 
    write_tsv(paste0("data/enrich/luma-id-comm-type-for-alluvial-", algrthm, ".tsv"))
 
  healthy_enr %>% left_join(healthy_comm_info %>% select(com_id, symbol), by = c("commun" = "com_id")) %>% 
    select(ID, symbol, community_type)%>% unique() %>% 
    write_tsv(paste0("data/enrich/healthy-id-comm-type-for-alluvial-", algrthm, ".tsv"))
  
  ### Enriched terms by comm
  luma_enr %>% group_by(commun, community_type) %>% tally() %>% 
    rename(community_id = commun, terms = n) %>%
    write_tsv(paste0("data/enrich/comm-enriched-terms-", algrthm, ".tsv"))
  
})
conds <- c("healthy", "luma")
all_enrich <- lapply(conds, function(cond) {
  all_enrich <- lapply(algorithms, function(algrthm) {
    alluvial_data <- read_tsv(paste0("data/enrich/",cond, "-id-comm-type-for-alluvial-", algrthm, ".tsv"))
    alluvial_data$algrthm <- algrthm
    alluvial_data$cond <- cond
    return(alluvial_data)
  })  
})


all_enrich <- bind_rows(all_enrich)
all_enrich %>% select(algrthm, symbol, community_type, cond) %>% distinct() %>% 
  group_by(community_type, algrthm, cond) %>% tally() %>% arrange(algrthm)

# community_type algrthm             cond        n
# <chr>          <chr>               <chr>   <int>
# 1 -              fast_greedy         healthy    14
# 2 -              fast_greedy         luma       20
# 3 *              fast_greedy         luma        9
# 4 -              infomap             healthy    47
# 5 -              infomap             luma       20
# 6 *              infomap             healthy     1
# 7 *              infomap             luma       16
# 8 -              leading_eigenvector healthy    18
# 9 -              leading_eigenvector luma       20
# 10 *              leading_eigenvector healthy     1
# 11 *              leading_eigenvector luma        9
# 12 -              multi_level         healthy    17
# 13 -              multi_level         luma       20
# 14 *              multi_level         luma        9