library(readr)
library(dplyr)
library(ggplot2)
library(parallel)
library(IRanges)

tads <- read_tsv("data/tads/GSE118712_MCF7_40kb_TADs.bed", col_names = c("chr", "start", "end"))
tads <- tads %>% mutate(size = end - start)
tads$chr <- sub("chr", "", tads$chr)

png(paste0("figures/tads/size-density.png"), width = 600, height = 400)
ggplot(tads, aes(x = size)) +
  scale_x_log10() +
  geom_density() + 
  theme_light()
dev.off()  

quantile(tads$size)
# 0%      25%      50%      75%     100% 
# 80000   440000   720000  1240000 19600000 

genes <- read_tsv("data/communities/luma-intra-vertices.tsv", col_types = cols_only(
  ensemblID = col_character(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  community = col_integer()
))

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                         col="black", sep=0.5, ...)
  {
    height <- 1
    if (is(xlim, "IntegerRanges"))
      xlim <- c(min(start(xlim)), max(end(xlim)))
      bins <- disjointBins(IRanges(start(x), end(x) + 1))
      
      plot.new()
    
      plot.window(xlim, c(0, max(bins)*(height + sep)))
    
      ybottom <- bins * (sep + height) - height
      
        rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
      
        title(main)
      
        axis(1)
}

chrs <- c(as.character(1:22), "X")

tads_comm_counts <- mclapply(X = chrs, 
  mc.cores = 75,
  FUN = function(ch){
    
    tads_chr <- tads %>% filter(chr == ch)
    genes_chr <- genes %>% filter(chr == ch)
    comms_chr <- unique(genes_chr$community)
    
    comms_info <- lapply(comms_chr, function(comm){
      gcomm <- genes_chr %>% filter(community == comm)
      return(list(comm = comm, start = min(gcomm$start), end =  max(gcomm$end), ngenes = nrow(gcomm)))
    })
    comms_info <- bind_rows(comms_info)
    comms_info <- comms_info %>% mutate(size = end - start)
    
    commranges <- IRanges(start = comms_info$start, end =  comms_info$end, names = comms_info$comm)
    plotRanges(commranges)
    tadranges <- IRanges(start = tads_chr$start, end =  tads_chr$end)
    plotRanges(tadranges)
    
    tadOverlaps <- as.data.frame(findOverlaps(commranges, tadranges, type = "within"))
       
    if(nrow(tadOverlaps) > 0) {
      tadOverlaps <- tadOverlaps %>% mutate(comm = comms_info[tadOverlaps$queryHits, "comm"] %>% unlist(use.names = F),
                                            init_space = comms_info[tadOverlaps$queryHits, "start"] %>% unlist(use.names = F) -
                                              tads_chr[tadOverlaps$subjectHits, "start"] %>% unlist(use.names = F),
                                            end_space =  tads_chr[tadOverlaps$subjectHits, "end"] %>% unlist(use.names = F) -
                                              comms_info[tadOverlaps$queryHits, "end"] %>% unlist(use.names = F), 
                                            comm_size = comms_info[tadOverlaps$queryHits, "size"] %>% unlist(use.names = F),
                                            tad_size =  tads_chr[tadOverlaps$subjectHits, "size"] %>% unlist(use.names = F),
                                            genes_comm = comms_info[tadOverlaps$queryHits, "ngenes"] %>% unlist(use.names = F))
      tadOverlaps$chr <- ch
    }
    return(tadOverlaps)
  })

tads_comm_counts <- bind_rows(tads_comm_counts)
sum(duplicated(tads_comm_counts$comm))
length(unique(genes$community))
write_tsv(tads_comm_counts, file = "data/tads/tads_comm_overlaps.tsv")
