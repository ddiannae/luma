library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

conds <- c("healthy", "luma")
label_conds <- c("Healthy", "Luma")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
colors <- c("#C7C7C7", "#FF9E4A")

all_comm_info <- lapply(conds, function(cond){
  comm_info <- read_tsv(paste0("data/", cond, "-communities-info.tsv"))
  comm_info$cond <- label_conds[cond]
  return(comm_info)
})
all_comm_info <- bind_rows(all_comm_info)
all_comm_info %>% group_by(cond) %>% filter( order > 200) 

p <- all_comm_info %>% 
  group_by(cond, gr = cut(order, 
                          breaks = c(0, 5, 50, 100, 200, 400, 800),
                          labels = c("1-5", "6-49", "50-99", "100-199", "200-399", "400-799"))) %>% 
  summarise(nnodes = n()) %>%
  arrange(as.numeric(gr)) %>% 
  ggplot( aes(y = nnodes, x = gr, fill=cond, color = cond)) +
  geom_bar(position=position_dodge(width=1), stat="identity") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_base()  +
  theme(legend.title = element_blank()) +
  #facet_wrap(~cond, nrow = 2) +
  xlab("Number of nodes in community") +
  ylab("Frequency")

png(paste0("figures/communities/community-order-by-cond.png"), width = 800, height = 800)
p
dev.off()  
