library(readr)
library(dplyr)
library(parallel)
library(IRanges)
library(ggplot2)

ctcfs <- read_tsv("data/ctcfs_in_intra_luma_50k.tsv", col_types = cols_only(
                  id = col_character(),
                  chr = col_character(),
                  start = col_integer(),
                  end = col_integer(),
                  type = col_character()))
ctcfs_hits <- read_tsv("data/ctcfs_hits_in_intra_luma.tsv")

genes <- read_tsv("data/luma-intra-vertices.tsv", col_types = cols_only(
  ensemblID = col_character(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  community = col_integer()
))
genes <- genes %>% dplyr::rename("id" = "ensemblID")
## get ranges for intergen ctcfs
intergen_ctcfs <- ctcfs %>% filter(type == "intergen")

r = 50000
## get the set of ctcfs associated to a community and its classification
## for that community
intra_comm <- unique(genes$community)
ctcfs_in_communities <- mclapply(X = intra_comm, 
         mc.cores = 75,
         FUN = function(comm){
           gene_comm <- genes %>% filter(community == comm)
           promoters <- gene_comm %>% mutate(end = start + 5000, start = start - 5000)
           extended_regions <- gene_comm %>% mutate(start = start - 5000, end = end)
           gene_bodies <- gene_comm %>% mutate(start = start + 5000) %>% filter(end > start)
           
           #ctcfs in genes or promoters
           gene_comm_ctcfs <- ctcfs %>% semi_join((ctcfs_hits %>% 
                                  semi_join(gene_comm, by = c("sequence" = "id"))), 
                               by = c("id" = "ctcf"))
           
           #intergen ctcfs in the community diameter and +-r
           s <- min(gene_comm$start)
           e <- max(gene_comm$end) 
           #ctcfs in the extended retion
           intergen_comm_ctcfs <- intergen_ctcfs %>% filter(chr == gene_comm$chr[1], 
                                      start >= s - r, end <= e + r)
           
           gctcfranges <- IRanges(start = gene_comm_ctcfs$start, end = gene_comm_ctcfs$end, 
                                  names = gene_comm_ctcfs$id)
           ictcfranges <- IRanges(start = intergen_comm_ctcfs$start, end = intergen_comm_ctcfs$end, 
                              names = intergen_comm_ctcfs$id)
           granges <- IRanges(start = gene_bodies$start, end = gene_bodies$end,
                              names = gene_bodies$id)
           pranges <- IRanges(start = promoters$start, end = promoters$end, names = promoters$id)
           eranges <- IRanges(start = extended_regions$start, end = extended_regions$end, 
                              names = extended_regions$id)
           
           # Classify gene ctcfs for this community
           gene_comm_ctcfs <- gene_comm_ctcfs %>% mutate(phits = countOverlaps(gctcfranges, pranges, type = "within"), 
                                             ghits = countOverlaps(gctcfranges, granges, type = "within"))

           gene_comm_ctcfs <- gene_comm_ctcfs %>% mutate(type = case_when(phits > 0 ~ "promoter",
                                                                  ghits > 0 ~ "gene"))
           # Distance to the extended region
           dist <- distanceToNearest(ictcfranges, eranges)
           intergen_comm_ctcfs$near_distance <- unlist(elementMetadata(dist))
           intergen_comm_ctcfs <- intergen_comm_ctcfs %>% 
             mutate(type = if_else(near_distance < r, "inter-less50", "inter-plus50"))
           
           comm_ctcfs  <- bind_rows(intergen_comm_ctcfs %>% select(id, chr, start, end, type),
                      gene_comm_ctcfs  %>% select(id, chr, start, end, type))

           return(comm_ctcfs)
      })
names(ctcfs_in_communities) <- intra_comm

## divide in slices and count
ctcfs_comm_counts <- mclapply(X = intra_comm, 
                                 mc.cores = 75,
                                 FUN = function(comm){
                                   
    comm_ctcfs <- ctcfs_in_communities[[as.character(comm)]]
    gene_comm <- genes %>% filter(community == comm)
    s <- min(gene_comm$start)
    e <- max(gene_comm$end) 
    
    #from gene start to gene end
    sls <- seq(s, e, by = r)
    
    #from gene start to gene end
    if(e - s > r) {
      #build all slices
      slices <- data.frame(start = c(s-r, sls, e),
                           end = c(s, sls[2:length(sls)], e, e+r))
    }else{
      slices <- data.frame(start = c(s-r, s, e),
                           end = c(s, e, e+r))
    }

    sranges <- IRanges(start = slices$start, end = slices$end)
    cranges <- IRanges(start = comm_ctcfs$start, end = comm_ctcfs$end, names = comm_ctcfs$id)
    #get ctcfs ovelaps 
    over_ctcfs <- as.data.frame(findOverlaps(cranges, sranges, type = "any")) %>% 
      mutate(id = comm_ctcfs$id[queryHits]) %>% select(-queryHits) %>% 
      dplyr::rename(slice = subjectHits) %>% left_join(comm_ctcfs, by = "id") %>% 
      mutate(slice_type = case_when(slice == 1 ~ "out_start",
                                    slice == length(sranges) ~ "out_end",
                                    TRUE ~ "mid"))
    if(nrow(over_ctcfs) > 0){
      over_ctcfs$comm = comm
    }
    return(list(overlap_ctcfs = over_ctcfs, 
            info_slices = data.frame(community = comm, 
                                      slices = length(sranges), 
                                      start = s,
                                      end = e,
                                      diameter = e - s,
                                      ctcfs_in_comm = nrow(comm_ctcfs),
                                      genes_in_comm = nrow(gene_comm)
                                     )))
})

info_slices <- plyr::ldply(lapply(ctcfs_comm_counts, "[[", "info_slices"))
ctcfs_comm_counts <-  plyr::ldply(lapply(ctcfs_comm_counts, "[[", "overlap_ctcfs"))

sum_by_types <- ctcfs_comm_counts %>% group_by(type, slice_type) %>% tally()

ggplot(sum_by_types, aes(x = slice_type,  y = n, fill = type)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() 

mid_ctcfs <- ctcfs_comm_counts %>% filter(slice_type == "mid") %>% 
  group_by(comm, type) %>% tally() %>%
  inner_join(info_slices, by = c("comm" = "community")) %>%
  mutate(ctcf_count = n/(slices-2), slice_type = "mid") %>% 
  select(comm, type, slice_type, ctcf_count)

outer_ctcfs <- ctcfs_comm_counts %>% filter(slice_type != "mid") %>% 
  group_by(comm, type, slice_type) %>% tally(name = "ctcf_count")

all_ctcfs <- bind_rows(mid_ctcfs, outer_ctcfs)
count_by_types <- all_ctcfs %>% group_by(type, slice_type) %>% 
  summarise(sum_ctcf = sum(ctcf_count))

ggplot(count_by_types, aes(x = slice_type,  y = sum_ctcf, fill = type)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() 

length(unique(ctcfs_comm_counts$comm))
# [1] 424

exp_assort <- read_tsv("data/assortativity/luma-exp-assortativity.tsv")
exp_assort <- exp_assort %>% dplyr::rename(comm = com, assort = totalfrac) %>% 
  select(comm, assort)

all_ctcfs_by_comm <- all_ctcfs %>% group_by(comm, slice_type) %>% 
  summarise(sum_ctcf = sum(ctcf_count)) %>% 
  tidyr::spread(key = slice_type, value = sum_ctcf, fill = 0) %>%
  mutate(type = case_when(out_end >= mid | out_start >= mid  ~ "out",
                          mid > out_end & mid > out_start ~ "mid",
                                                      TRUE ~ "NA")) %>%
  inner_join(exp_assort, by = "comm") %>% 
  inner_join(info_slices %>% select(-start, -end), by = c("comm" = "community")) 

ggplot(all_ctcfs_by_comm, aes(x = diameter, color = type))+
  geom_density() +
  scale_x_log10() +
  theme_bw() +
  scale_color_viridis_d() 

ggplot(all_ctcfs_by_comm, aes(x = assort, color = type))+
  geom_density() +
  theme_bw() +
  scale_color_viridis_d() 

ggplot(all_ctcfs_by_comm, aes(x = ctcfs_in_comm, color = type))+
  geom_density() +
  theme_bw() +
  scale_color_viridis_d() 

ggplot(all_ctcfs_by_comm, aes(x = genes_in_comm, color = type))+
  geom_histogram(bins = 50) +
  theme_bw() +
  scale_color_viridis_d() 

ggplot(all_ctcfs_by_comm, aes(x = type, fill = type)) +
  geom_bar() +
  theme_bw() +
  scale_fill_viridis_d() 

 