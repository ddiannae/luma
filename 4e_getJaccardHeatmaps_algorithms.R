library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ComplexHeatmap)
library(circlize)
library(janitor)

conds <- c("healthy", "luma")
label_conds <- c("Healthy", "LumA")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
names(colors) <- conds
algorithms <- c("fast_greedy", "infomap" ,"leading_eigenvector", "multi_level")
labels_alg <- c("Fast Greedy", "Infomap", "Leading Eigenvector", "Louvain")
names(labels_alg) <- algorithms

col_fun = colorRamp2(c(0, 0.1+1e-6,  0.5 , 1), 
                     c("white","#0D0887FF", "#CC4678FF",  "#F0F921FF"))

biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt") %>% clean_names() %>%
  select(gene_stable_id, hgnc_symbol)

alg_pairs <- as.data.frame(combn(algorithms, 2), stringsAsFactors = F) 

l <- lapply(conds, function(cond){
  m <- apply(alg_pairs, MARGIN = 2, FUN=function(ap)  {
      alg1 <- ap[1]
      alg2 <- ap[2]
      
      jaccs <- read_tsv(paste0("data/communities/", cond, "-", alg1, "-", alg2, ".tsv"))
      total1 <- max(jaccs$comm1)
      total2 <- max(jaccs$comm2)
      non_zero <- jaccs %>% filter(jaccard != 0)
      label_size <- c("(0.0 - 0.05)", "[0.05, 0.10)", "[0.10, 0.25)", 
                      "[0.25, 0.50)", "[0.50, 0.75)", "[0.75, 1.0)", "1.0")
      
      p <- non_zero %>%
        mutate(gr = case_when(
          .$jaccard < 0.05 ~ label_size[1],
          .$jaccard < 0.10 ~ label_size[2],
          .$jaccard < 0.25 ~ label_size[3],
          .$jaccard < 0.50 ~ label_size[4],
          .$jaccard < 0.75 ~ label_size[5],
          .$jaccard < 1.0 ~ label_size[6],
          TRUE ~ label_size[7]
        )) %>% 
        group_by(gr) %>%
        summarise(nnodes = n()) %>%
        mutate(gr = factor(gr, levels = label_size)) %>%
        arrange(gr) %>% 
        ggplot( aes(y = nnodes, x = gr)) +
        geom_bar(aes(color = cond, fill = cond), 
                 position=position_dodge(width=1), stat="identity") +
        geom_text(aes(label = nnodes), vjust = -0.2, size = 3)  +
        scale_fill_manual(values = colors[cond]) +
        scale_color_manual(values = colors[cond]) +
        theme_base()  +
        xlab("Jaccard Index") +
        ylab("Frequency") +
        scale_y_continuous(expand = expansion(add = 100)) +
        theme(legend.position = "none", plot.background=element_blank(), 
              plot.title = element_text(size = 12), 
              text = element_text(size=10)) + 
        ggtitle(paste0(labels_alg[alg1], " (", total1, ") vs ", 
                       labels_alg[alg2], " (", total2, ")" ), 
                subtitle = "Non zero values")
      
      png(paste0("figures/communities/", cond, "-", alg1, "-", alg2,
                 "-bars.png"), height = 3,
          width = 5, units = "in", res = 300)
      print(p)
      dev.off()
      
      one_to_one <- jaccs %>% filter(jaccard == 1) 
      
      comm1_info  <- read_tsv(paste0("data/communities/", cond, "-communities-info-", 
                                      alg1, ".tsv"))
      comm2_info  <- read_tsv(paste0("data/communities/", cond, "-communities-info-", 
                                     alg2, ".tsv"))
      
      jaccard_matrix <- jaccs %>% 
        anti_join(one_to_one, by="comm1") %>% 
        anti_join(one_to_one, by="comm2")  %>% 
        left_join(comm1_info %>% select(com_id, pg_gene, order),
                  by=c("comm1" ="com_id")) %>%
        left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
        select(-pg_gene) %>% rename(gene_comm1 = hgnc_symbol, order_comm1 = order) %>%
        left_join(comm2_info %>% select(com_id, pg_gene, order),
                  by=c("comm2" ="com_id")) %>%
        left_join(biomart, by =c("pg_gene"= "gene_stable_id")) %>% 
        select(-pg_gene) %>% rename(gene_comm2 = hgnc_symbol,  order_comm2 = order) %>%
        arrange(comm1, comm2) %>%
        mutate(c1 = paste0(gene_comm1, " (", order_comm1, ")"),
               c2 = paste0(gene_comm2, " (", order_comm2, ")")) %>%
        pivot_wider(id_cols = c1, names_from = c2, values_from = jaccard)
      
      comm1_names <- jaccard_matrix %>% select(c1) %>% unlist(use.names = F)
      jaccard_matrix <- jaccard_matrix %>% select(-c1) %>% as.matrix()
      rownames(jaccard_matrix) <- comm1_names
      
      total1 <- nrow(jaccard_matrix)
      total2 <- ncol(jaccard_matrix)
      
      # if(total2 > total1) {
      #   jaccard_matrix <- t(jaccard_matrix)
      #   tmp_alg <- alg1
      #   alg1 <- alg2
      #   alg2 <- tmp_alg
      #   tmp_total <- total1
      #   total1 <- total2
      #   total2 <- tmp_total
      # }
      
      ht <- Heatmap(jaccard_matrix, name = "Jaccard\n index", col = col_fun, 
                    cluster_rows = F, cluster_columns = F,
                    column_title = paste0( labels_alg[alg2],  " (", total2, ")"),
                    row_title = paste0( labels_alg[alg1],  " (", total1, ")"),
                    column_title_gp = gpar(fontsize = 18),
                    row_names_gp = gpar(fontsize = 8), 
                    column_names_gp = gpar(fontsize = 8),
                    row_title_gp = gpar(fontsize = 18),
                    heatmap_legend_param=list(title_gp=gpar(fontsize=6), 
                                              labels_gp=gpar(fontsize=4)))
      
      
      png(paste0("figures/communities/",  cond, "-", alg1, 
                 "-", alg2, "-heatmap.png"), height = 5,
          width = 6, units = "in", res = 300)
        draw(ht)
      dev.off()  
  })
})
