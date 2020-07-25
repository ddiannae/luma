library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

all_lesions <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv", 
                        col_types = cols_only(id = col_integer(), type = col_character(), 
                                              q = col_double())) 
min_q <- min(all_lesions$q)
max_q <- max(all_lesions$q)

all_genes_in_lesions <- read_tsv("data/luma-gistic/all-genes-in-lesions.tsv")

comm_enrich <- read_tsv("data/enrich-universe/comm-enriched-terms.tsv",
                        col_types = cols_only(community_id = col_double(), terms = col_double()))
### Solo para comunidades con más de 30 terms
comm_enrich <- comm_enrich %>% filter(terms > 30)
comm_membership <- read_tsv("data/communities/luma-communities.tsv") %>% 
  rename(ensembl_id = ensemblID, community_id = community)
comm_membership <- comm_membership %>% semi_join(comm_enrich, by = "community_id")

comm_membership <- comm_membership %>% left_join(all_genes_in_lesions, by = "ensembl_id") %>%
  left_join(read_tsv("data/network-tables/luma-20127-vertices.tsv",
                     col_types = cols_only(ensemblID = col_character(), symbol = col_character(), 
                                           chr = col_character(), band = col_character())),
            by = c("ensembl_id" = "ensemblID")) %>% 
  mutate(chr = factor(chr, levels = as.character(c(seq(1:22), "X"))))

comm_membership <- comm_membership %>% left_join(all_lesions, by = c("type", "id")) %>%
  mutate(lq = -log10(q), tlq = ifelse(type == "del", lq*-1, lq))

#### Color palettes
col_fun = colorRamp2(c(log10(min_q), 0, -log10(min_q)),
                     c("blue", "white", "red"))
col_fun = colorRamp2(c(-8, -5, 0, 5, 8),
                     c("#0000AA", "#0000AA", "white", "#AA0000", "#AA0000"))

chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
                     "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
                     "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
                     "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9")

names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                            "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21")

### from the inter lfc enrichment plot
comms <- c(316, 158, 146, 197, 183, 230, 36, 575, 91)

## falta la 2,3
## 11 heatmaps for each community
com_name <- comms[3]
comm <- comm_membership %>% filter(community_id == com_name) %>% arrange(chr) 
comm_matrix <- comm %>%  select(tlq) %>%  as.matrix() 
rownames(comm_matrix) <- comm %>% select(symbol) %>% unlist()


chrs <- as.data.frame(comm %>% select(chr) %>% 
                        rename(Chr = chr))


ha <- rowAnnotation(df = chrs, 
                    name = "Chr", show_annotation_name = F,
                    col = list(Chr = chromosomes.pal[as.character(chrs$Chr)]), 
                    width = unit(0.5, "cm"),
                    annotation_legend_param = list(
                      title = "Chr", title_gp = gpar(fontsize = 14), grid_height = unit(0.5, "cm")))

ht <- Heatmap(comm_matrix, cluster_rows = F,  col = col_fun, na_col = "white",
              heatmap_legend_param = list( title_gp = gpar(fontsize = 14),
                title = "-log10(q)", at = c(-8, -5, 0, 5, 8), legend_height = unit(4, "cm"),
                labels = c("", "5-Del", "0", "5-Amp", "")
              ),
              name = "-log10(q)", rect_gp = gpar(col = "gray", lwd = 2),
              row_names_side = "left", width = unit(1, "cm"),
              right_annotation = ha,
              show_column_names = F, show_row_names = T)

draw(ht[101:103,], heatmap_legend_side = "right", annotation_legend_side = "right")

