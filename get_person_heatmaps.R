library(readr)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

cond <- "luma"
pearson_matrix <- read_csv("data/pearson/luma_pearson.csv", 
                         col_types = cols(
                           src = col_character(),
                           dst = col_character(),
                           src_gstart = col_double(),
                           dst_gstart = col_double(),
                           src_chrom = col_character(),
                           dst_chrom = col_character(),
                           pearson = col_double()))

pearson_matrix <- pearson_matrix %>% 
  filter(src_chrom == "1")

pearson_matrix <- pearson_matrix %>% 
  pivot_wider(id_cols = src, names_from = dst,values_from = pearson)
srcs <- pearson_matrix[, 1] %>% unlist(use.names = F)

pearson_matrix <- pearson_matrix %>% select(-src) %>% mutate(newcol = NA) %>% 
  select(newcol, everything()) 
colnames(pearson_matrix)[1] <- srcs[1]

pearson_matrix <- pearson_matrix %>% as.matrix()
rownames(pearson_matrix) <- srcs
#pearson_matrix <- pearson_matrix + t(pearson_matrix)
diag(pearson_matrix) <- 1

col_pal <- colorRamp2(c(-1, -0.5,  0, 0.5, 1), 
                      c("blue4", "blue3", "white","red3", "red4"))
ht <- Heatmap(pearson_matrix, col = col_pal, na_col = "white", 
              cluster_rows = F,  name = "Pearson coeff",
              cluster_columns = F, 
              show_column_names = F, show_row_names = F)

png(paste0("figures/", cond, "-pearson-heatmap.png"),
    units="in", width=20, height=20, res=75)
print(ht)
dev.off()
