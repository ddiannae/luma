library(readr)
library(dplyr)

annot <- read.delim("data/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col.names = c("ensemblID", "chr", "start", "end", "band"), stringsAsFactors = F)
annot <- annot[!duplicated(annot$ensemblID), ]

conds <- c("luma", "healthy")
n <- lapply(conds, function(cond){
  cat("Reading network ", cond, "\n")
  net <- read_tsv(paste0("data/", cond, ".1e8.sif"), col_names = c("source", "MI", "target"))
  
  cat("Is not unsorted: ", !is.unsorted(rev(net$MI)), "\n")
  net <- net %>% arrange(desc(MI))
  net$cond <- cond
  net$row_num <- 1:nrow(net)
  
  colnames(annot) <-  c("source", "source_chr", "source_start", "source_end", "source_band")
  net <- merge(annot, net)
  colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_band")
  net <- merge(annot, net)
  net <- net %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                        interaction_type = case_when(inter == F & source_band == target_band ~ "Intra-Cytoband",
                                                     inter == F & source_band != target_band ~ "Inter-Cytoband",
                                                     TRUE ~ "Trans"),
                        distance = if_else(inter == F, pmax(source_start, target_start) - 
                                             pmin(source_start, target_start), as.integer(NaN)))
  
  annot.symbol <- read.delim("data/Biomart_EnsemblG94_GRCh38_p12.txt",
                             colClasses = c("character",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "character" ), 
                             stringsAsFactors = F, 
                             col.names = c("ensemblID",  "NULL",  "NULL",  "NULL", "NULL", "NULL",  "symbol"))
  
  annot.symbol <- annot.symbol[!duplicated(annot.symbol$ensemblID), ]
  targets <- net %>% select(target, target_chr, target_start, target_end, target_band)
  sources <- net %>% select(source, source_chr, source_start, source_end, source_band)
  colnames(targets) <- c("ensemblID", "chr", "start", "end", "band")
  colnames(sources)  <- c("ensemblID", "chr", "start", "end", "band")
  vertices <- rbind(targets, sources)
  vertices <- vertices[!duplicated(vertices$ensemblID), ]
  vertices <- vertices %>% inner_join(annot.symbol, by = "ensemblID")
  
  housekeeping <- read.delim("data/housekeeping.txt", header = F, 
                             colClasses = c("character", "NULL"))
  housekeeping <- housekeeping[, 1, drop= T]
  housekeeping <- stringr::str_trim(housekeeping)
  housekeeping <- tolower(housekeeping)
  vertices <- vertices %>% mutate(isHK = if_else(tolower(symbol) %in% housekeeping, T, F))
  
  tfs <- read.delim("data/tfs_2018.txt", header = F,
                    col.names = c("ensemblID"))
  vertices <- vertices %>% mutate(isTF = if_else(ensemblID %in% tfs$ensemblID, T, F))
  
  # kariotype <- read.delim("/media/ddisk/transpipeline-data/annotations/chromosome.band.hg38.txt",
  #                         header = T, stringsAsFactors = F)
  # names(kariotype)[1] <- "chr"
  # tbands <- kariotype %>% mutate(arm = unlist(lapply(as.vector(strsplit(name, split = "")), "[[", 1)))
  # tbands <- tbands %>% group_by(chr, arm) %>% summarise(t = max(name))
  # tbands <- tbands %>% group_by(chr) %>% summarise(t1 = min(t), t2 = max(t))
  # tbands$chr <- stringr::str_replace(tbands$chr, "chr", "")
  # vertices <- vertices %>% inner_join(tbands, by = "chr")
  # vertices$isExtreme <- vertices$band == vertices$t1 | vertices$band == vertices$t2
  # vertices$isExtreme <- as.integer(vertices$isExtreme)
  vertices$isHK <- as.integer(vertices$isHK)
  vertices$isTF <- as.integer(vertices$isTF)
  
  net <- net %>% 
    select(source, target, MI, distance, row_num, interaction_type) %>%
    arrange(desc(MI)) 
  
  # vertices <- vertices %>% select(-t1, -t2)
  
  write.table(vertices, file = paste0("data/network-tables/", cond, "-vertices.tsv") , 
              quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(net, file = paste0("data/network-tables/", cond, "-interactions.tsv"),
              quote = F, row.names = F, col.names = T, sep = "\t")
})

