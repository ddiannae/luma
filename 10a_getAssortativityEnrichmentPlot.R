library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

luma_exp <- read_tsv("data/assortativity/luma-exp-assortativity.tsv")

luma_chr <- read_tsv("data/assortativity/luma-chr-assortativity.tsv",
                     col_types = cols_only(community_id = col_double(),
                                           totalfrac = col_double()))
luma_chr <- luma_chr %>% filter(totalfrac < 1)

luma_assort <- luma_chr %>% inner_join(luma_exp, by = "community_id", suffix = c("_chr", "_exp"))

comm_info <- read_tsv("data/communities/luma-communities-info.tsv", 
                      col_types = cols_only(com_id = col_double(), order = col_double(), 
                                            pg_gene = col_character())) %>%
    left_join(read_tsv("data/network-tables/luma-20127-vertices.tsv",
                     col_types = cols_only(ensemblID = col_character(), symbol = col_character())),
            by = c("pg_gene" = "ensemblID"))

comm_enrich <- read_tsv("data/enrich-universe/comm-enriched-terms.tsv",
                        col_types = cols_only(community_id = col_double(), terms = col_double()))

luma_plot <- luma_assort %>% inner_join(comm_info, by = c("community_id" = "com_id")) %>%
  left_join(comm_enrich, by = "community_id") %>%  mutate(terms = if_else(is.na(terms), 0, terms))
luma_plot[luma_plot$community_id == 158, "symbol"] <- "FOXM1"

p <- ggplot(luma_plot, aes(x = totalfrac_chr, y = mean_diff_exp)) + 
  geom_point(aes(size = order, color = terms)) +
  geom_label(aes(label = ifelse(terms > 30, as.character(symbol), NA)),
            colour = "darkslategrey", size = 6,
            nudge_y = -0.20, nudge_x = 0.0, label.padding = unit(0.15, "lines"),) +
  geom_hline(yintercept = 0, linetype="dashed", color = "gray") +
  theme_base() +
  labs(size = "Nodes in community",
       color = "Enriched terms") +
  xlab("Fraction of intra-chromosomal links") +
  ylab("Mean LFC") +
  ylim(c(-3, 2.0)) +
  scale_color_gradient(low="#CCFFCC", high="#003333") +
  scale_size_continuous(range = c(1, 15), breaks = c(1, 5, 20, 50, 100, 150, 200)) +
  theme(text = element_text(size = 20), axis.title = element_text(size = 25))

png(paste0("figures/communities/intercomms-lfc-enrichment.png"), width = 1200, height = 600)
p
dev.off()
