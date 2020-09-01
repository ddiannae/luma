library(readr)
library(dplyr)

luma_enr <- read_tsv("data/enrich-universe/luma-go-filtered-enrichments.tsv", 
                     col_types = cols_only(ID = col_character(), Description = col_character(),
                                           p.adjust = col_double(), commun = col_integer()))

luma_info_comm <- read_tsv("data/communities/luma-communities-info.tsv") %>% 
  left_join(read_tsv("data/communities/luma-intra-communities.tsv"), by = c("com_id" = "community_id")) %>%
  mutate(intra = !is.na(diameter))

luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(ensemblID = col_character(), symbol = col_character(),
                                                isTF = col_logical(), community = col_integer()))
all_tfs <- luma_vertices %>% semi_join(luma_enr %>% select(commun) %>% distinct(), 
                            by = c("community" = "commun")) %>%  filter(isTF == TRUE) %>% 
  select(symbol, community) %>% inner_join(luma_info_comm, by = c("community" = "com_id")) %>%
  left_join(luma_vertices %>% select(ensemblID, symbol), by = c("pg_gene" = "ensemblID"), 
            suffix = c("_tf", "_pg"),)
  
comm_tfs <- all_tfs %>% group_by(community) %>% summarise(ntfs = n(), order = first(order), 
                                              intra = first(intra), pg = first(symbol_pg)) %>% 
  mutate(dif = order-ntfs) %>% arrange(desc(dif))

#     community  ntfs order intra pg       dif
# <dbl> <int> <dbl> <lgl> <chr>  <dbl>
#  1       169     7   177 FALSE NAA30    170
#  2       281     5   110 FALSE UTP15    105
#  3       146     2   103 FALSE RPL35    101
#  4       158     4    80 FALSE NUSAP1    76
#  5       183     3    26 FALSE CCL5      23
#  6       197     1    23 FALSE MX1       22
#  7       349     1    23 FALSE ESAM      22
#  8       316     1    20 FALSE COL5A2    19
#  9       230     8    20 FALSE NR4A1     12
# 10        57     1     7 FALSE IGLL5      6
# 11        36     1     6 FALSE TYROBP     5
# 12        14     4     6 FALSE PER3       2
# 13       134     6     7 TRUE  HOXB3      1
# 14         4     6     6 TRUE  HOXA5      0
# 15       475     5     5 TRUE  HOXC6      0

lapply(comm_tfs %>% select(community) %>% unlist(), function(community_id) {
  all_tfs %>% filter(community == community_id) %>% select(symbol_tf) %>%
    write_tsv(paste0("data/tfs/tfs-community-", community_id, ".txt"), col_names = FALSE)
})


