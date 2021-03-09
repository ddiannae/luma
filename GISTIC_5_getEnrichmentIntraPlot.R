library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggthemes)

all_enrichments <- read_tsv("data/luma-gistic/community-enrichment.tsv")
gistic_peaks <- read_tsv("data/luma-gistic/all-lesions-conf-99.tsv") 
gistic_peaks <- gistic_peaks %>% mutate(id = paste0(id, "_", type))

vertices <- read_tsv("data/network-tables/luma-20127-vertices.tsv")

algrthm <- "multi_level"
comm_info <- read_tsv(paste0("data/communities/luma-communities-info-",algrthm, ".tsv"), 
                      col_types = cols_only(com_id = col_double(), order = col_double(), 
                                            pg_gene = col_character())) %>%
  left_join(vertices %>% select(ensemblID, symbol),
            by = c("pg_gene" = "ensemblID"))

comm_expr <- read_tsv(paste0("data/assortativity/luma-exp-assortativity-", algrthm, ".tsv")) %>% 
  select(community_id, diffraction, mean_diff_exp)

comm_chr <- read_tsv(paste0("data/assortativity/luma-chr-assortativity-", algrthm, ".tsv"),
                     col_types = cols_only(community_id = col_double(),
                                           diffraction = col_double()))
gene_comm <- read_tsv(paste0("data/communities/luma-communities-", algrthm, ".tsv"))

comms <- unique(gene_comm$community)

luma_assort <- comm_chr %>% 
  inner_join(comm_expr, by = "community_id", suffix = c("_chr", "_exp"))

comm_info <- comm_info %>% inner_join(luma_assort, by = c( "com_id" = "community_id"))

intra_comm <- comm_info %>% filter(diffraction_chr == 1)

vertices <- vertices %>% inner_join(gene_comm, by = "ensemblID")

comm_starts <- lapply(intra_comm$com_id, function(comm){
  vertices %>% filter(community == comm) %>%  slice(which.min(start)) %>% 
    select(community, start, chr)
})
comm_starts <- bind_rows(comm_starts)

comm_info <- comm_info %>% inner_join(comm_starts, by = c("com_id" = "community"))
comm_info %>% arrange(desc(chr))
comm_info$chrom <- as.numeric(ifelse(comm_info$chr == "X", 23, comm_info$chr))
comm_info <- comm_info %>% arrange(chrom, start)
comm_info$xpos <- 1:nrow(comm_info)
comm_info <- comm_info %>% left_join(all_enrichments, by = c("com_id" = "comm")) 

comm_info <- comm_info %>% rename(lesion_id = ID)
comm_info <- comm_info %>% left_join(gistic_peaks, by =c("lesion_id" = "id"), 
                                     suffix = c("_comm", "_peak"))

comm_info$chr_comm <- factor(comm_info$chr_comm, levels = as.character(c(1:22, "X")))

## only the visible ones
chr_labs <- c(as.character(1:12), "", as.character(14:17), "", as.character(19:20), "", "22", "X")
names(chr_labs) <-  as.character(c(1:22, "X"))

p <- ggplot(comm_info, aes(x = xpos, y = mean_diff_exp)) + 
  geom_bar(stat = "identity", aes(color = type_log10_q)) +
  geom_point(aes(size = order),  fill = NA, shape=19) +
  theme_base() +
  labs(size = "Nodes in \ncommunity", color = "Peak\n-log10(q)") +
  xlab("") +
  ylab("Mean LFC") +
  scale_color_gradient2(midpoint=0, low="turquoise", mid="white",
                        high="maroon1", space ="Lab",  na.value = "grey80",
                        limits = c(-2, 2), breaks = c(-2, 0, 2)) +
  scale_size_continuous(range = c(0, 3), breaks = c(1, 2, 5, 50)) +
  theme(text = element_text(size = 12), axis.title = element_text(size = 12),
        plot.background=element_blank()) +
  facet_grid(~chr_comm, scales = 'free_x', switch = "x", space = 'free_x',
             labeller = labeller(chr_comm = chr_labs)) + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 8), 
        strip.background = element_rect(fill = "grey90", colour = "white"),
        axis.line.x = element_blank()) +
  ggtitle("Intra-chromosomal communities with Amp-Del peaks enrichment")

png(paste0("figures/gistic/gistic-enrichments.png"), units="in",
    width=8, height=4, res=300)
print(p)
dev.off()
