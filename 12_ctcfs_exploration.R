library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)

annot <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col_names = c("ensembl_id", "chr", "start", "end", "karyoband"),
                    col_types = cols(col_character(), col_character(), col_integer(), 
                                     col_integer(), col_character()),
                    skip = 1)
genes <- read_tsv("data/genes_in_exp_matrix.txt", col_names = "ensembl_id")

annot <- annot %>% semi_join(genes)

ctcfs <- read_tsv("data/non_gene_overlaping_ctcfbs.tsv", 
                  col_types = cols(col_character(), col_integer(), col_integer(),
                                   col_character(), col_double()))
ctcfs <- ctcfs[!duplicated(ctcfs), ]

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ctcfs_counts <- ctcfs %>% group_by(chr) %>% tally()
ctcfs_counts$chr <- factor(ctcfs_counts$chr, levels = as.character(c(seq(1:22), "X")))

png(filename = "figures/ctcf_count.png", width = 1200, height = 600)
ggplot(ctcfs_counts, aes(x = chr, y = n, fill = chr)) +
  geom_bar(stat = "identity") +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_manual(values = getPalette(23))
dev.off()

chrs <- as.character(c(seq(1:22), "X"))

gb_ctfs <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  gs_in_c <- annot %>% filter(chr == ch) %>% 
    select(-karyoband) %>%
    dplyr::rename(id = ensembl_id) %>%
    mutate(type = "Gene")
  cs_in_c <- ctcfs %>% filter(chr == ch) %>%
    select(-score)%>%
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

ceros <- gb_ctfs %>% filter(genes == 0) %>% group_by(chr) %>% tally()
ceros$empty = "Empty"
noceros <- gb_ctfs %>% filter(genes != 0) %>% group_by(chr) %>% tally()
noceros$empty = "Non Empty"

ctcfs_counts <- bind_rows(ceros, noceros)

png(filename = "figures/ctcf_count_empty.png", width = 1200, height = 600)
ggplot(ctcfs_counts, aes(x = chr, y = n)) +
  geom_bar(
    aes(color = empty, fill = empty),
    stat = "identity", position = position_stack()
  ) +
  theme_few(base_size = 20) +
  theme(legend.title = element_blank()) +
  ylab("CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_tableau() +
  scale_color_tableau() 
dev.off()

noceros <- gb_ctfs %>% filter(genes != 0) 
png(filename = "figures/ctcf_genes_between_nonempty.png", width = 1200, height = 600)
ggplot(noceros, aes(x = chr, y = genes, fill = chr)) +
  geom_boxplot() +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("Genes between CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_manual(values = getPalette(23))
dev.off()

png(filename = "figures/ctcf_distance_all.png", width = 1200, height = 600)
ggplot(gb_ctfs, aes(x = chr, y = distance, fill = chr)) +
  geom_boxplot() +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("Distance between CTCF binding sites") +
  xlab("Chromosome") +
  ylim(0, 5e5) +
  scale_fill_manual(values = getPalette(23))
dev.off()

write_tsv(gb_ctfs, path = "data/all_ctcfs.tsv")