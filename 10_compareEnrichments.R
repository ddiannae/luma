library(readr)
library(dplyr)

healthy_enr <- read_tsv("./data/enrich-universe/healthy-go-enrichments.tsv")
luma_enr <- read_tsv("./data/enrich-universe/luma-go-enrichments.tsv")

#### No filter for size
healthy_enr <- healthy_enr %>% filter(p.adjust < 0.005)
luma_enr <- luma_enr %>% filter(p.adjust < 0.005)

healthy_enr <- healthy_enr %>% select(ID, Description) %>% unique()
luma_enr <- luma_enr %>% select(ID, Description) %>% unique()

shared <- healthy_enr %>% semi_join(luma_enr)
luma_only <- luma_enr %>% anti_join(healthy_enr)
healthy_only <- healthy_enr %>% anti_join(luma_enr)

write_tsv(luma_only %>% select(ID), "data/enrich-universe/luma-only-GO.tsv")
write_tsv(healthy_only %>% select(ID), "data/enrich-universe/healthy-only-GO.tsv")
write_tsv(shared %>% select(ID), "data/enrich-universe/shared-only-GO.tsv")

healthy_kegg <- read_tsv("./data/enrich-universe/healthy-kegg-enrichments.tsv")
luma_kegg <- read_tsv("./data/enrich-universe/luma-kegg-enrichments.tsv")

healthy_kegg <- healthy_kegg %>% filter(p.adjust < 0.005)
luma_kegg <- luma_kegg %>% filter(p.adjust < 0.005)

healthy_kegg <- healthy_kegg %>% select(ID, Description) %>% unique()
luma_kegg <- luma_kegg %>% select(ID, Description) %>% unique()

shared <- healthy_kegg %>% semi_join(luma_kegg)
luma_only <- luma_kegg %>% anti_join(healthy_kegg)
healthy_only <- healthy_kegg %>% anti_join(luma_kegg)

write_tsv(luma_only %>% select(ID), "data/enrich-universe/luma-only-KEGG.tsv")
write_tsv(healthy_only %>% select(ID), "data/enrich-universe/healthy-only-KEGG.tsv")
write_tsv(shared %>% select(ID), "data/enrich-universe/shared-only-KEGG.tsv")
