library(readr)
library(stringr)
library(dplyr)
library(IRanges)

annot <- read_tsv("data/Biomart_EnsemblG94_GRCh38_p12.txt", skip = 1,
                  col_names = c("ensembl_id", "chr", "start", "end", "GC", "type", "symbol"))

amps <- read_tsv("data/luma-gistic/amp_genes.conf_99.txt",
                 col_names = FALSE)

amps_info <- amps[1:4, ] %>% t() 
rownames(amps_info) <- NULL
colnames(amps_info) <- c("cytoband", "q", "residual_q", "wide_peak")
amps_info <- amps_info[c(-1, -46),]
amps_info <- amps_info %>% as_tibble() %>% 
  mutate(q =  as.numeric(q), residual_q = as.numeric(residual_q), 
         id = rownames(.), neg_log10_q = -log10(q))

dels <- read_tsv("data/luma-gistic/del_genes.conf_99.txt",
                 col_names = FALSE)

dels_info <- dels[1:4, ] %>% t() 
rownames(dels_info) <- NULL
colnames(dels_info) <- c("cytoband", "q", "residual_q", "wide_peak")
dels_info <- dels_info[c(-1, -37),]
dels_info <- dels_info %>% as_tibble() %>% 
  mutate(q =  as.numeric(q), residual_q = as.numeric(residual_q), 
         id = rownames(.),  neg_log10_q = -log10(q))

amps_info <- amps_info %>% mutate(pos = str_split(wide_peak, ":"), 
                                  chr = str_remove(unlist(lapply(pos, "[[", 1)), "chr"),
                                  se = str_split(unlist(lapply(pos, "[[", 2)), "-"),
                                  start = as.integer(unlist(lapply(se, "[[", 1))),
                                  end = as.integer(unlist(lapply(se, "[[", 2)))) %>% 
  select(-pos, -se) %>% mutate(type = "amp")

dels_info <- dels_info %>% mutate(pos = str_split(wide_peak, ":"), 
                                  chr = str_remove(unlist(lapply(pos, "[[", 1)), "chr"),
                                  se = str_split(unlist(lapply(pos, "[[", 2)), "-"),
                                  start = as.integer(unlist(lapply(se, "[[", 1))),
                                  end = as.integer(unlist(lapply(se, "[[", 2)))) %>% 
  select(-pos, -se) %>% mutate(type = "del")

chrs <- as.character(c(seq(1:22), "X"))

all_genes_lesions <- parallel::mclapply(X = chrs, mc.cores = 23, FUN = function(ch){
  chr_amps <- amps_info %>% filter(chr == ch)
  chr_genes <- annot %>% filter(chr == ch)
  chr_dels <- dels_info %>% filter(chr == ch)
  
  aranges <- IRanges(start = chr_amps$start, end = chr_amps$end, names = chr_amps$id)
  dranges <- IRanges(start = chr_dels$start, end = chr_dels$end, names = chr_dels$id)
  granges <- IRanges(start = chr_genes$start, end = chr_genes$end, names = chr_genes$ensembl_id)

  genes <- as.data.frame(findOverlaps(granges, aranges, type = "within"))
  
  amps <- genes %>% 
    mutate(ensembl_id = chr_genes$ensembl_id[queryHits],
           id = chr_amps$id[subjectHits], type = "amp") %>% select(ensembl_id, id, type)
  
  genes <- as.data.frame(findOverlaps(granges, dranges, type = "within"))
  
  dels <- genes %>% 
    mutate(ensembl_id = chr_genes$ensembl_id[queryHits],
           id = chr_dels$id[subjectHits], type = "del") %>% select(ensembl_id, id, type)
  
  return(bind_rows(amps, dels))
})

bind_rows(all_genes_lesions)   %>% write_tsv(path = "data/luma-gistic/all-genes-in-lesions.tsv")
bind_rows(amps_info, dels_info) %>% select(id, type, everything()) %>%
  write_tsv(path = "data/luma-gistic/all-lesions-conf-99.tsv")
