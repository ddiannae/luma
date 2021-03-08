library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(janitor)

algorithms <- c( "fast_greedy", "infomap", "leading_eigenvector","multi_level")

m <- lapply(algorithms, function(algrthm) {
  
  luma_exp <- read_tsv(paste0("data/assortativity/luma-exp-assortativity-", algrthm,".tsv")) %>% 
    select(community_id, diffraction, mean_diff_exp)
  
  luma_chr <- read_tsv(paste0("data/assortativity/luma-chr-assortativity-",algrthm ,".tsv"),
                       col_types = cols_only(community_id = col_double(),
                                             diffraction = col_double()))
  luma_chr <- luma_chr %>% filter(diffraction < 1)
  
  luma_assort <- luma_chr %>% inner_join(luma_exp, by = "community_id", suffix = c("_chr", "_exp"))
  
  
  comm_info <- read_tsv(paste0("data/communities/luma-communities-info-", algrthm, ".tsv"), 
                        col_types = cols_only(com_id = col_double(), order = col_double(), 
                                              pg_gene = col_character())) %>%
      left_join(read_tsv("data/network-tables/luma-20127-vertices.tsv",
                       col_types = cols_only(ensemblID = col_character(), symbol = col_character())),
              by = c("pg_gene" = "ensemblID"))
  
  comm_enrich <- read_tsv(paste0("data/enrich/luma-go-enrichments-", algrthm, ".tsv")) %>% 
    clean_names() %>% select(id, commun) %>% group_by(commun) %>% tally() %>% 
    rename("nterms" = "n") 
  
  luma_plot <- luma_assort %>% inner_join(comm_info, by = c("community_id" = "com_id"))%>%
    left_join(comm_enrich, by = c("community_id" = "commun")) 
  luma_plot[is.na(luma_plot$nterms), "nterms"] <- 0
  
  luma_plot %>% filter(nterms > 20) %>% 
    dplyr::select(community_id) %>% write_tsv(paste0("data/enrich/communities-more-20-terms-",algrthm, ".tsv"))
  
  luma_plot %>% select(community_id, symbol, order, diffraction_chr, diffraction_exp, mean_diff_exp, nterms) %>%
    rename(size = order, name = symbol, chr_assortativity = diffraction_chr, 
           exp_assprtativity = diffraction_exp, enriched_terms = nterms) %>% 
    write_tsv(paste0("data/enrich/luma-assort-terms-info-",algrthm, ".tsv"))
  
  p <- ggplot(luma_plot, aes(x = diffraction_chr, y = diffraction_exp)) + 
    geom_point(aes(size = order, color = nterms)) +
    geom_text(aes(label = ifelse(nterms >20, as.character(symbol), NA)),
              colour = "black", size = 1.3, check_overlap = FALSE, fontface = "bold") +
    #geom_text(aes(label = ifelse(nterms <20 & nterms >0, as.character(symbol), NA)),
    #          colour = "black", size = 1, check_overlap = FALSE) +
    geom_hline(yintercept = 0.0, linetype="dashed", color = "gray") +
    geom_vline(xintercept = 0.0, linetype="dashed", color = "gray") +
    theme_base() +
    labs(size = "Nodes in \ncommunity",
         color = "Enriched \nterms") +
    xlab("Chromosomal assortativity") +
    ylab("Expression assortativity") +
    scale_color_gradient(low="lightblue1", high="steelblue4") +
    scale_size_continuous(range = c(1, 7), breaks = c(1, 5, 20, 50, 100, 150, 200)) +
    theme(text = element_text(size = 6), 
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 4),
          plot.background=element_blank(), legend.key.size = unit(0.6, "line"),
          legend.spacing.y = unit(0.01, 'cm'))
  
  png(paste0("figures/communities/intercomms-lfc-enrichment-",  algrthm, ".png"), 
      units="in", width=4, height=2.5, res=300)
  print(p)
  dev.off()

})
