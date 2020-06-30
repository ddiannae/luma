library(readr)
library(dplyr)

ctcfs <- read_tsv("data/all_ctcfs.tsv", col_types = cols_only(
              id = col_character(),
              chr = col_character(),
              start = col_double(),
              end = col_double()
            ))

luma_vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv", 
                          col_types = cols_only(
                              ensemblID = col_character(),
                              chr = col_character(),
                              start = col_double(),
                              end = col_double()
                          ))
luma_vertices <- luma_vertices %>% rename(id = ensemblID)

chrs <- as.character(c(seq(1:22), "X"))

gb_ctfs <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  gs_in_c <- luma_vertices %>% filter(chr == ch) %>% 
    mutate(type = "Gene")
  cs_in_c <- ctcfs %>% filter(chr == ch) %>%
    mutate(type = "CTCF")
  
  all <- bind_rows(gs_in_c, cs_in_c) %>% arrange(start) %>%
    mutate(gene_sum = cumsum(type == "Gene"))
  cs_in_c <- all %>% filter(type == "CTCF")
  
  ## Check how many genes are between two ctcf binding sites
  ngenes <- lapply(seq(1:nrow(cs_in_c)), function(i){
    cs <- unname(unlist(cs_in_c[i, "gene_sum"]))
    start <- unname(unlist(cs_in_c[i, "start"]))
    pre_cs <- unname(unlist(cs_in_c[i-1, "gene_sum"]))
    pre_end <- unname(unlist(cs_in_c[i-1, "end"]))
    
    return(list(genes = cs - pre_cs, distance = start - pre_end))
  })
  cs_in_c$genes <- c(cs_in_c$gene_sum[1], unlist(lapply(ngenes, "[[", "genes")))
  cs_in_c$distance <- c(cs_in_c$start[1], unlist(lapply(ngenes, "[[", "distance")))
  
  return(cs_in_c)
})

gb_ctfs <- bind_rows(gb_ctfs) %>% select(-type)
extreme_ctcfs <- gb_ctfs %>% group_by(chr, gene_sum) %>% filter(
  end == min(end) | end == max(end)
)
write_tsv(extreme_ctcfs, path = "data/extreme_ctcfs_in_luma.tsv")
