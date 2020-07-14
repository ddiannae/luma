library(readr)
library(dplyr)
library(parallel)
library(IRanges)

ctcfs <- read_tsv("data/ctcfs_in_intra_luma_50k.tsv", col_types = cols_only(
                  id = col_character(),
                  chr = col_character(),
                  start = col_integer(),
                  end = col_integer()))

genes <- read_tsv("data/luma-intra-vertices.tsv", col_types = cols_only(
  ensemblID = col_character(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  community = col_integer()
))

chrs <- as.character(c(seq(1:22), "X"))
ctcfs_ranges <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  chr_ctcfs <- ctcfs %>% filter(chr == ch)
  return(IRanges(start = chr_ctcfs$start, end = chr_ctcfs$end, names = chr_ctcfs$id))
})
names(ctcfs_ranges) <- chrs

w = 10000
r = 50000
ctcfs_in_communities <- mclapply(X = unique(genes$community), 
         mc.cores = 75,
         FUN = function(comm){
           gene_comm <- genes %>% filter(community == comm)
           #gene start
           s <- min(gene_comm$start)
           #gene end
           e <- max(gene_comm$end) 
           #from region start to gene start
           ss <- seq(s-r , s, by = w)
           #from gene end to region end
           es <- seq(e, e+r, by = w)
           #from gene start to gene end
           if(e - s > w) {
             mes <- seq(s, e, by = w)
             #build all slices
             slices <- data.frame(start = c(ss[1:length(ss) -1], mes, es[1:length(es) -1]),
                                  end = c(ss[2:length(ss)], mes[2:length(mes)], e, es[2:length(es)]))
           }else{
             slices <- data.frame(start = c(ss[1:length(ss) -1], s, es[1:length(es) -1]),
                                  end = c(ss[2:length(ss)], e, es[2:length(es)]))
           }
           
           slices <- IRanges(start = slices$start, end = slices$end)
           #get ctcfs ovelaps (any)
           ctcfs_comm <- as.data.frame(findOverlaps(ctcfs_ranges[[gene_comm$chr[1]]], 
                                                    slices))
           if(nrow(ctcfs_comm) > 0){
              ctcfs_comm$comm <- comm
           }
           return(list(ctcfs_comm = ctcfs_comm, 
                       total_slices = data.frame(community = comm, 
                                                 slices = length(slices), 
                                                 diameter = e - s)))
         })

comm_slices <- plyr::ldply(lapply(ctcfs_in_communities, "[[", "total_slices"))
ctcfs_in_communities <-  plyr::ldply(lapply(ctcfs_in_communities, "[[", "ctcfs_comm"))

length(unique(ctcfs_in_communities$comm))
# [1] 438

outss <- r/w
ctcfs_in_communities <- ctcfs_in_communities %>% group_by(comm, subjectHits) %>% tally()
ctcfs_in_communities <- ctcfs_in_communities %>%
  mutate(total_slices = comm_slices %>% filter(community == comm) %>% 
          select(slices) %>% unlist(use.names= FALSE), 
         type = case_when(subjectHits <= outss ~ "start",
                          subjectHits >= total_slices - outss ~ "end",
                                    TRUE ~ "mid")) %>%
  select(-total_slices)

ctcfs_counts <- ctcfs_in_communities %>% group_by(comm, type) %>% 
  summarise(ctcfs = sum(n)) %>% tidyr::spread("type", "ctcfs", fill = 0)
