library(readr)
library(dplyr)
library(limma)
library(Glimma)
library(edgeR)

conds <- c("healthy", "luma")

all_matrices <- lapply(conds, function(cond){
  matrix <- read_tsv(paste0("data/", cond, "_cpm10_arsyn.tsv"))
  return(list(matrix = matrix[,-1], nrows = nrow(matrix), ncols = ncol(matrix[,-1]), 
              genes = matrix[,1], cond = cond))
})

matrices <- lapply(all_matrices, "[[", 1)
nrows <- unlist(lapply(all_matrices, "[[", 2))
ncols <- unlist(lapply(all_matrices, "[[", 3))
genes <- unlist(all_matrices[[1]][4], use.names = FALSE)
conds <- unlist(lapply(all_matrices, "[[", 5))

M <- cpm(data.matrix(bind_cols(matrices)), log = T)
rownames(M) <- genes
targets <- data.frame(id = colnames(M), group = factor(unlist(mapply(rep, conds, ncols))))

## MUST BE TRUE ##
all(nrows == length(genes))

## DISEÑO ##
dmatrix <- model.matrix(~1 + targets$group , data = targets)

## AJUSTE ##
fit <- lmFit(M, dmatrix)
head(fit$coefficients)

fit <- eBayes(fit)

fit$fdr <- apply(fit$p.value, 2, p.adjust, method="fdr")
p <- 1-fit$fdr 
fit$B <- log(p/(1-p))

topTable(fit, coef=ncol(dmatrix), adjust="BH")

#Glimma plot
dt <- decideTests(fit, lfc = 2) 
summary(dt)
glMDPlot(path = "data", fit, counts = M, groups = targets$group, status = dt)

genesFULL <- bind_cols(
  ensembl_id = rownames(fit$coefficients),
  lfc = fit$coefficients[, 2],
  avg_log2_cpm_exp = fit$coefficients[, 1] + fit$coefficients[, 2],
  p_value = fit$p.value[, 2],
  FDR = fit$fdr[, 2],
  B = fit$B[, 2]
)
genesFULL$avg_exp <- rowMeans(matrices[[2]])
genesFULL$avg_log2_exp <- rowMeans(log2(matrices[[2]]))
genesFULL$avg_cpm_exp <- rowMeans(cpm(matrices[[2]], log = FALSE))

genesFULL <- genesFULL %>% mutate(pass = abs(lfc) >= 0.5 | p_value < 0.05, 
                                  lfc = ifelse(pass, lfc, 0),
                                  p_value = ifelse(pass, p_value, 0),
                                  FDR = ifelse(pass, FDR, 0),
                                  B = ifelse(pass, B, 0)) %>% select(-pass)

write_tsv(genesFULL, paste0("data/", conds[2], "-deg-ebayes.tsv"))

n <- 20127  
vertices <- read_tsv(file = paste0("data/network-tables/", conds[2], "-", n, "-vertices.tsv"))
vertices <- vertices %>%  inner_join(genesFULL, by = c("ensemblID" = "ensembl_id"))

write.table(vertices, file = paste0("data/network-tables/", conds[2], "-", n, "-vertices.tsv") , 
            quote = F, row.names = F, col.names = T, sep = "\t")
