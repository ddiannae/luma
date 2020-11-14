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

library(RColorBrewer)
library(viridis)
bluecols <- brewer.pal(9, 'YlGnBu')
newcol <- colorRampPalette(bluecols)

pie(rep(1, ncols), col = bluecols2, border = NA, labels = NA)


comm_info <- read_tsv("data/communities/luma-communities-info.tsv")
plotRanges <- function(x, heights = heights, cols = colors, xlim=x, main=deparse(substitute(x)),
                         sep=0.5, ...)
  {
    height <- 1
    ncols <- max(cols)
    newcolors <- inferno(ncols)[cols]
    legend_image <- as.raster(matrix(rev(inferno(ncols)), ncol=1))
    if (is(xlim, "IntegerRanges"))
      xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    heights <- heights*3
    heights <- ifelse(heights>1, 1, heights)
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + heights, col=newcolors, ...)
    title(main)
    axis(1)
    plot.new()
    text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
    rasterImage(legend_image, 0, 0, 1,1)
      
}

chrs <- c(as.character(1:22), "X")

tads_comm_counts <- mclapply(X = chrs, 
  mc.cores = 75,
  FUN = function(ch){
    
    tads_chr <- tads %>% filter(chr == ch)
    genes_chr <- genes %>% filter(chr == ch)
    # comms_chr <- unique(genes_chr$community)
    # 
    # comms_info <- lapply(comms_chr, function(comm){
    #   gcomm <- genes_chr %>% filter(community == comm)
    #   return(list(comm = comm, start = min(gcomm$start), end =  max(gcomm$end), ngenes = nrow(gcomm)))
    # })
    # comms_info <- bind_rows(comms_info)
    # comms_info <- comms_info %>% mutate(size = end - start)
    
    generanges <- IRanges(start = genes_chr$start, end =  genes_chr$end, names = genes_chr$ensemblID)
    #plotRanges(generanges)
    tadranges <- IRanges(start = tads_chr$start, end =  tads_chr$end)
    #plotRanges(tadranges)
    
    tadOverlaps <- as.data.frame(findOverlaps(generanges, tadranges, type = "any"))
       
    if(nrow(tadOverlaps) > 0) {
      # tadOverlaps <- tadOverlaps %>% mutate(comm = comms_info[tadOverlaps$queryHits, "comm"] %>% unlist(use.names = F),
      #                                       init_space = comms_info[tadOverlaps$queryHits, "start"] %>% unlist(use.names = F) -
      #                                         tads_chr[tadOverlaps$subjectHits, "start"] %>% unlist(use.names = F),
      #                                       end_space =  tads_chr[tadOverlaps$subjectHits, "end"] %>% unlist(use.names = F) -
      #                                         comms_info[tadOverlaps$queryHits, "end"] %>% unlist(use.names = F), 
      #                                       comm_size = comms_info[tadOverlaps$queryHits, "size"] %>% unlist(use.names = F),
      #                                       tad_size =  tads_chr[tadOverlaps$subjectHits, "size"] %>% unlist(use.names = F),
      #                                       genes_comm = comms_info[tadOverlaps$queryHits, "ngenes"] %>% unlist(use.names = F))
      tadOverlaps <- tadOverlaps %>% mutate(ensemblID =  genes_chr[tadOverlaps$queryHits, "ensemblID"] %>% unlist(use.names = F))
    }
    return(tadOverlaps)
  })
tads_comm_counts <- bind_rows(tads_comm_counts)
tads_comm_counts <- tads_comm_counts %>% left_join(genes, by = "ensemblID") %>% group_by(community) %>% 
  summarise(chr = min(chr), genes_in_tads = n_distinct(ensemblID), span_tads = n_distinct(subjectHits)) 
tads_comm_counts <- tads_comm_counts %>% left_join(comm_info %>% select(-chr), by = c("community" = "com_id"))
tads_comm_counts <- tads_comm_counts %>% mutate(gene_dif = order - genes_in_tads)
tads_comm_counts %>% arrange(span_tads)

comm_ranges_info <- mclapply(X = chrs, 
                     mc.cores = 75,
                     FUN = function(ch){
                       
                       genes_chr <- genes %>% filter(chr == ch)
                       comms_chr <- unique(genes_chr$community)

                       comms_info <- lapply(comms_chr, function(comm){
                         gcomm <- genes_chr %>% filter(community == comm)
                         return(list(comm = comm, start = min(gcomm$start), end =  max(gcomm$end)))
                       })
                       comms_info <- bind_rows(comms_info)
                       comms_info <- comms_info %>% mutate(comm_diam = end - start)
                       return(comms_info)
})
comm_ranges_info <- bind_rows(comm_ranges_info)

comm_ranges_info <- comm_ranges_info %>% left_join(tads_comm_counts, by = c("comm"="community"))


crom_1 <- comm_ranges_info%>% filter(chr == "1")
crom_1$ngenes <- crom_1$order/max(crom_1$order)
commranges = IRanges(start = crom_1$start, end = crom_1$end)
plotRanges(commranges, crom_1$ngenes, cols = crom_1$span_tads)

showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}

cols <- rev(inferno(70))
showCols(cl= cols, bg="gray33", rot=30, cex=0.75)

layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(1:20, 1:20, pch = 19, cex=2, col = inferno(20))

tads_chr1 = tads %>% filter(chr == "1")
tadranges = IRanges(tads_chr1$start, tads_chr1$end)
plotOInlyRanges(tadranges)


plotOInlyRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col="black", ...)
  title(main)
  axis(1)
  
}


plotRanges(tadranges, cols = c(1))
legend_image <- as.raster(matrix(inferno(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
layout(1, width = 1,height = 1)
plot.new()


write_tsv(tads_comm_counts, file = "data/tads/tads_gene_count_overlap.tsv")
