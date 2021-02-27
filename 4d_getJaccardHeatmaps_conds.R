library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ComplexHeatmap)
library(circlize)

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

m <- lapply(algorithms, function(alg)  {
   
    jaccs <- read_tsv(paste0("data/communities/healthy-luma-", alg, ".tsv"))
    total1 <- max(jaccs$healty)
    total2 <- max(jaccs$luma)
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
      geom_bar( position=position_dodge(width=1), stat="identity") +
      geom_text(aes(label = nnodes), vjust = -0.2)  +
    #  scale_fill_manual(values = colors) +
    #  scale_color_manual(values = colors) +
      theme_base()  +
      xlab("Jaccard Index") +
      ylab("Frequency") +
      theme(legend.position = "none", plot.background=element_blank()) + 
      ggtitle(paste0(label_conds[1], " (", total1, ") vs ", 
                     label_conds[2], " (", total2, ")" ), 
              subtitle = "Non zero values")
    
    png(paste0("figures/communities/jaccard-index-", cond, "-", alg1, "-", alg2,
               "-bars.png"),  width = 600, height = 350)
    print(p)
    dev.off()
    
    one_to_one <- jaccs %>% filter(jaccard == 1) 
    
    jaccard_matrix <- jaccs %>% 
      anti_join(one_to_one, by="comm1") %>% 
      anti_join(one_to_one, by="comm2")  %>% 
      arrange(comm1, comm2) %>%
      pivot_wider(id_cols = comm1, names_from = comm2, values_from = jaccard)
    
    comm1_names <- jaccard_matrix %>% select(comm1) %>% unlist(use.names = F)
    jaccard_matrix <- jaccard_matrix %>% select(-comm1) %>% as.matrix()
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
                  row_title_gp = gpar(fontsize = 18),
                  heatmap_legend_param=list(title_gp=gpar(fontsize=12), 
                                            label_gp=gpar(fontsize=10)))
    
    
    png(paste0("figures/communities/jaccard-index-",  cond, "-", alg1, 
               "-", alg2, "-heatmap.png"), height = 8,
        width = 8, units = "in", res = 100)
    draw(ht)
    dev.off()  
  })
})