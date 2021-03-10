library(dplyr)
library(readr)
library(janitor)
library(ggplot2)
library(ggthemes)

path <- "/datos/ot/diana/subtipos_mama_2018/networks/"
biomart <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt") %>% clean_names()
biomart <- biomart %>% 
  select(gene_stable_id, chromosome_scaffold_name, gene_start_bp, hgnc_symbol)

inters <- c(10, 100, 1000)

conds <- c("healthy","luma")
label_conds <- c("Healthy", "LumA")
names(label_conds) <- conds
colors <- c("#e3a098", "#a32e27")
names(colors) <- conds


dfs <- lapply(conds, function(cond) {
  ll <- lapply(inters, function(inter) {
    net <- read_tsv(paste0(path, cond, "_", inter, "k.tsv"), 
                    col_names = c("from", "MI", "to")) 
    
    cnames <- colnames(biomart)
    colnames(biomart) <- paste0("from_",colnames(biomart))
    net <- net %>% inner_join(biomart, by = c("from" = "from_gene_stable_id"))
    
    colnames(biomart) <- cnames
    colnames(biomart) <- paste0("to_",colnames(biomart))
    
    net <- net %>% inner_join(biomart, by = c("to" = "to_gene_stable_id"))
    net <- net %>% mutate(type = if_else(from_chromosome_scaffold_name == to_chromosome_scaffold_name, "cis", "trans"))
    
    sets <- seq(1, nrow(net), inter*100)
    
    ll <- lapply(sets, function(start) {
      end <- start + inter*100-1
      subnet <- net[start:end, ]
      cis <- sum(subnet$type == "cis", na.rm = T)
      return(cis/(sum(subnet$type == "trans", na.rm = T)+cis))
    })
    ll <- unlist(ll, use.names = F)
    df <- data.frame(start = sets, frac = ll, end = sets + (inter*100-1))
    df$cond <- cond
    return(df)
  })
})

dfs_df <- bind_rows(dfs)
dfs_df <- dfs_df %>% filter((start != 1 & end != 10000) |
                 (start != 1 & end != 100000)) %>%
          mutate(cond_label = label_conds[cond])


p <- ggplot(dfs_df, aes(y = frac, x = end, color = cond)) +
  geom_point(size = 2) +
  scale_x_log10(name="Number of Interactions", ) +
  theme_base() +
  theme(  text = element_text(size=8)) +
  ylab("Fraction of -cis interactions") +
  scale_color_manual(name="Condition", labels = label_conds, values = colors) 
  
png(paste0("figures/cis-networks-interactions.png"), 
    units="in", width=6, height=2.5, res=300)
print(p)
dev.off()
