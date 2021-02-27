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

col_fun = colorRamp2(c(0, 0.1+1e-6,  0.5 , 1), 
                     c("white","#0D0887FF", "#CC4678FF",  "#F0F921FF"))

all_jaccards <- lapply(algorithms, function(alg)  {
    jaccs <- read_tsv(paste0("data/communities/healthy-luma-", alg, ".tsv"))
    jaccs$alg <- alg
    total_healthy <- max(jaccs$healthy)
    total_luma <- max(jaccs$luma)
    jaccs <- jaccs %>% filter(jaccard != 0)
    return(list(nonzero = jaccs, total_luma = total_luma,
                total_healthy = total_healthy, alg = alg))
})

non_zero <- bind_rows(lapply(all_jaccards, "[[", "nonzero"))
labels <-  bind_rows(lapply(all_jaccards, "[",
                            c("total_luma", "total_healthy", "alg")))
labels <- labels %>% mutate(
  label = paste0(labels_alg[alg], ". Healthy (", total_healthy, 
                 ") - Luma (", total_luma, ")")
)

labels_alg <- labels %>% select(label) %>% unlist(use.names = F)
names(labels_alg) <- labels$alg


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
  group_by(gr, alg) %>%
  summarise(nnodes = n()) %>%
  mutate(gr = factor(gr, levels = label_size)) %>%
  arrange(gr) %>% 
  ggplot( aes(y = nnodes, x = gr)) +
  geom_bar(position=position_dodge(width=1), stat="identity", 
           fill = "steelblue4") +
  geom_text(aes(label = nnodes), vjust = -0.2)  +
  theme_base()  +
  xlab("Jaccard Index") +
  ylab("Frequency") +
  scale_y_continuous(expand = expansion(add = 150)) +
  theme(legend.position = "none", plot.background=element_blank(),
        strip.text = element_text(size = 20)) + 
  facet_wrap(~alg, labeller = as_labeller(labels_alg)) 

png("figures/communities/jaccard-index-all-conds-all-algs-bars.png",
      width = 1200, height = 800)
print(p)
dev.off() 


# Which communities have jaccard index = 1?
one_to_one <- non_zero %>% filter(jaccard == 1) 
    
one_to_one %>% group_by(healthy_n, luma_n) %>% tally()
# healthy_n luma_n     n
# <dbl>  <dbl> <int>
#   1         2      2    98
#   2         3      3     1
