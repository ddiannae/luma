library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

conds <- c("healthy", "luma")
label_conds <- c("Healthy", "LumA")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
names(colors) <- conds
algorithms <- c("Fast Greedy", "Infomap", "Leading Eigenvector", 
                "Louvain")
names(algorithms) <- c( "fast_greedy", "infomap", "leading_eigenvector", "multi_level")
m <- lapply(conds, function(cond){
  print(cond)
  
  all_comm_info <- lapply(names(algorithms), function(algrthm){
    comm_info <- read_tsv(paste0("data/communities/", cond, "-communities-info-",
                                 algrthm, ".tsv"), col_types = cols(chr = col_character()))
    comm_info$algrthm <- algorithms[algrthm]
    comm_info$cond <- cond
    return(comm_info)
  })
  
  all_comm_info <- bind_rows(all_comm_info)
  label_size <- c("1-2","3-4", "5-9", "10-19", "20-49", "50-99", "100-199",
                  "200-500", ">500")
    
  p <- all_comm_info %>% 
    group_by(algrthm, gr = cut(order, 
                            breaks = c(0, 3, 5, 10, 20, 50, 100, 200, 500, 
                                       max(all_comm_info$order)+1),
                            labels = label_size,
                            right = FALSE)) %>% 
    summarise(nnodes = n()) %>%
    pivot_wider(names_from = "gr", values_from = "nnodes", values_fill = NA) %>%
    pivot_longer(cols = !algrthm, names_to = "gr", values_to = "nnodes") %>%
    mutate(gr = factor(gr, levels = label_size)) %>%
    arrange(gr) %>% 
    ggplot( aes(y = nnodes, x = gr, fill=cond, color = cond)) +
    geom_bar(position=position_dodge(width=1), stat="identity") +
    scale_fill_manual(values = colors[cond]) +
    scale_color_manual(values = colors[cond]) +
    theme_base()  +
    theme(legend.position = "none") +
    xlab("Number of nodes in community") +
    ylab("Frequency") +
    facet_wrap(~algrthm, ncol = 1, labeller = label_value)
  
  png(paste0("figures/communities/community-order-by-algorithm-",cond,"-bars.png"), 
      width = 600, height = 800)
  print(p)
  dev.off()  
  
  p <- all_comm_info %>% 
    filter(order >= 3) %>%
    ggplot( aes(y = order, x = algrthm, fill=cond, color = cond)) +
    geom_violin() +
    scale_fill_manual(values = colors[cond]) +
    scale_color_manual(values = colors[cond]) +
    theme_base()  +
    theme(legend.position = "none") +
    ylab("Number of nodes in community") 
  
  png(paste0("figures/communities/community-order-by-algorithm-", cond, "-violin.png"), 
      width = 600, height = 400)
  print(p)
  dev.off()  
  
})
