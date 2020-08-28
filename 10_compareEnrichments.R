library(readr)
library(dplyr)

#### Get exclusive enrichments
healthy_enr <- read_tsv("./data/enrich-universe/healthy-go-enrichments.tsv")
luma_enr <- read_tsv("./data/enrich-universe/luma-go-enrichments.tsv")

#### No filter for size
healthy_enr <- healthy_enr %>% filter(p.adjust < 0.005)
luma_enr <- luma_enr %>% filter(p.adjust < 0.005)

write_tsv(luma_enr, "data/enrich-universe/luma-go-filtered-enrichments.tsv")
write_tsv(healthy_enr, "data/enrich-universe/healthy-go-filtered-enrichments.tsv")

healthy_enr <- healthy_enr %>% select(ID, Description) %>% unique()
luma_enr <- luma_enr %>% select(ID, Description) %>% unique()

shared <- healthy_enr %>% semi_join(luma_enr)
luma_only <- luma_enr %>% anti_join(healthy_enr)
healthy_only <- healthy_enr %>% anti_join(luma_enr)

write_tsv(luma_only %>% select(ID), "data/enrich-universe/luma-only-GO.tsv")
write_tsv(healthy_only %>% select(ID), "data/enrich-universe/healthy-only-GO.tsv")
write_tsv(shared %>% select(ID), "data/enrich-universe/shared-only-GO.tsv")

healthy_kegg <- read_tsv("./data/enrich-universe/healthy-kegg-enrichments.tsv")
luma_kegg <- read_tsv("./data/enrich-universe/luma-kegg-enrichments.tsv")

healthy_kegg <- healthy_kegg %>% filter(p.adjust < 0.005)
luma_kegg <- luma_kegg %>% filter(p.adjust < 0.005)

healthy_kegg <- healthy_kegg %>% select(ID, Description) %>% unique()
luma_kegg <- luma_kegg %>% select(ID, Description) %>% unique()

shared <- healthy_kegg %>% semi_join(luma_kegg)
luma_only <- luma_kegg %>% anti_join(healthy_kegg)
healthy_only <- healthy_kegg %>% anti_join(luma_kegg)

write_tsv(luma_only %>% select(ID), "data/enrich-universe/luma-only-KEGG.tsv")
write_tsv(healthy_only %>% select(ID), "data/enrich-universe/healthy-only-KEGG.tsv")
write_tsv(shared %>% select(ID), "data/enrich-universe/shared-only-KEGG.tsv")

### Data for alluvial
luma_enr <- read_tsv("data/enrich-universe/luma-go-filtered-enrichments.tsv")
luma_intra_comm <- read_tsv("data/communities/luma-intra-communities.tsv")

luma_enr <- luma_enr %>% mutate(community_type = 
                                  ifelse(commun %in% luma_intra_comm$community_id, 
                                         "*", "-"))
intra_enrich <- luma_enr %>% filter(community_type == "*")
## 136 términos en 9 comunidades intra
inter_enrich <-  luma_enr %>% filter(community_type == "-")
## 792 términos en 20 comunidades inter

luma_comm_info <- read_tsv("data/communities/luma-communities-info.tsv")
luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(ensemblID = col_character(), symbol = col_character()))

luma_comm_info <- luma_comm_info %>% left_join(luma_vertices, by = c("pg_gene" = "ensemblID"))

luma_enr %>% left_join(luma_comm_info %>% select(com_id, symbol), by = c("commun" = "com_id")) %>% 
  select(ID, symbol, community_type)%>% unique() %>% write_tsv("data/enrich-universe/id-comm-type-for-alluvial.tsv")

### Enriched terms by comm
luma_enr %>% group_by(commun, community_type) %>% tally() %>% 
  rename(community_id = commun, terms = n) %>%
  write_tsv("data/enrich-universe/comm-enriched-terms.tsv")

minterms <- luma_enr %>% left_join(luma_comm_info %>% select(com_id, symbol), by = c("commun" = "com_id")) %>% 
  select(ID, symbol, Description, community_type, commun, p.adjust)%>% 
  group_by(commun) %>% slice_min(p.adjust)  %>% write_tsv("data/enrich-universe/id-comm-type-min-terms.tsv")
