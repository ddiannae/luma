library(readr)
library(dplyr)
library(parallel)
library(IRanges)
library(ggplot2)
library(dgof)

ctcfs <- read_tsv("data/ctcfs_in_intra_luma_50k_1000_100.tsv", col_types = cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  id = col_character(),
  score = col_double(),
  size = col_double(),
  type = col_character(),
  near_region = col_character(),
  near_distance = col_double()
))

ctcfs_hits <- read_tsv("data/ctcfs_hits_in_intra_luma_1000_100.tsv")

genes <- read_tsv("data/luma-intra-vertices.tsv", col_types = cols_only(
  ensemblID = col_character(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  community = col_integer()
))

genes <- genes %>% dplyr::rename("id" = "ensemblID")
## get ranges for intergen ctcfs
intra_comm <- unique(genes$community)
r = 50000

## divide in slices and count
## Slices of lenght r
ctcfs_comm_counts <- mclapply(X = intra_comm, 
                                 mc.cores = 75,
                                 FUN = function(comm){
                                   
    gene_comm <- genes %>% filter(community == comm)
    comm_ctcfs <- ctcfs %>% filter(near_region %in% gene_comm$id)
    
    ## to get the promoter region associated
    s <- min(gene_comm$start) - 1000
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
    # get ctcfs ovelaps. 
    # Slice one, the start of the community. The last slice is the end of the community. 
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

ggplot(sum_by_types, aes(x = slice_type,  y = n, fill = type, width = 0.8)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() 

### Average number of CTCFs per slice 
mid_ctcfs <- ctcfs_comm_counts %>% filter(slice_type == "mid") %>% 
  group_by(comm, type) %>% tally() %>%
  inner_join(info_slices, by = c("comm" = "community")) %>%
  mutate(ctcf_count = n/(slices-2), slice_type = "mid") %>% 
  select(comm, type, slice_type, ctcf_count)

outer_ctcfs <- ctcfs_comm_counts %>% filter(slice_type != "mid") %>% 
  group_by(comm, type, slice_type) %>% tally(name = "ctcf_count")

all_ctcfs <- bind_rows(mid_ctcfs, outer_ctcfs)
count_by_types <- all_ctcfs %>% dplyr::group_by(type, slice_type, comm) %>% 
  dplyr::summarise(sum_ctcf = sum(ctcf_count)) 

ggplot(count_by_types, aes(x = slice_type,  y = sum_ctcf, fill = type)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_fill_viridis_d() +
  xlab("") +
  ylab("Total number of CTCF bs")

length(unique(ctcfs_comm_counts$comm))
# [1] 404
# 416 con 1000 100

## Fill the missing values with zeroes
all_ctcfs_by_comm <- all_ctcfs %>% dplyr::group_by(comm, slice_type) %>%
  dplyr::summarise(sum_ctcf = sum(ctcf_count)) %>%
  tidyr::spread(key = slice_type, value = sum_ctcf, fill = 0) %>%
  tidyr::gather(key = "slice_type", value = "sum_ctcf", mid, out_end, out_start) %>%
  arrange(comm)

ggplot(all_ctcfs_by_comm, aes(y = sum_ctcf, x = slice_type)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_viridis_d()

ggplot(all_ctcfs_by_comm, aes(y = sum_ctcf, x = slice_type, color = slice_type))+
  geom_violin() +
  theme_bw() +
  scale_color_viridis_d()

png(paste0("figures/ctcfs/density-count-by-type-1000-100.png"), width = 800, height = 400)
ggplot(all_ctcfs_by_comm, aes(x = sum_ctcf, color = slice_type))+
  geom_density(show.legend=FALSE) +
  stat_density(aes(x=sum_ctcf, colour=slice_type), 
               geom="line",position="identity") +
  theme_bw(base_size = 20) +
  scale_color_viridis_d(name = "Position", labels = c("Medium", "End", "Start")) +
  xlab("Average CTCFs binding sites") +
  ylab("Density") 
dev.off()


mids <- all_ctcfs_by_comm %>% filter(slice_type == "mid") %>% ungroup() %>% select(sum_ctcf) %>% 
  unlist(use.names = F)
pmids <- ecdf(mids)
vals = seq(-3, 3, by=0.01) 
out_end <- all_ctcfs_by_comm %>% filter(slice_type == "out_end") %>% ungroup() %>% select(sum_ctcf) %>% 
  unlist(use.names = F)
pouts <- ecdf(out_end)
out_start <- all_ctcfs_by_comm %>% filter(slice_type == "out_start") %>% ungroup() %>% select(sum_ctcf) %>% 
  unlist(use.names = F)
pstart <- ecdf(out_start)

ks.test(out_end, mids, exact = F, simulate.p.value = T)


### count number of communities with greater outs
comms_gather <- all_ctcfs_by_comm %>% 
  tidyr::spread(key = slice_type, value = sum_ctcf, fill = 0)
comms_gather <- comms_gather %>% mutate(type = case_when((out_end > mid | out_start > mid) ~ "out",
                                                         (mid > out_end | mid > out_start ) ~ "mid",
                                                         TRUE ~ "equal")
                                        )
comms_gather %>% group_by(type) %>% tally()

# type      n
# 1 equal     1
# 2 mid     219
# 3 out     196