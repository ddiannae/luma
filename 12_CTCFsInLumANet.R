library(readr)
library(dplyr)

ctcfs <- read_tsv("data/all_ctcfs.tsv", col_types = cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  id = col_character(),
  score = col_double(),
  phits = col_double(),
  ghits = col_double(),
  ohits = col_double(),
  promoter = col_character(),
  gen = col_character(),
  other = col_character()))

luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(
                            ensemblID = col_character(),
                            chr = col_character()
                          ))
luma_vertices <- luma_vertices %>% rename(id = ensemblID)

ctcfs <- ctcfs %>% filter(gen %in% luma_vertices$id | 
                            promoter %in% luma_vertices$id)

write_tsv(ctcfs, path = "data/ctcfs_in_luma.tsv")
