###################################################################################
### Author: QuanZhang
### time: 2023-05-05
### Downstream analysis of miRNA
###################################################################################

## load library
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(tibble)
library(org.Mm.eg.db)
library(FactoMineR)
library(factoextra)
library(export)
library(AnnotationDbi)
library(readxl)

#### load data
mirna <- read.delim("./data/miRNAexpression/Readcount_TPM.xls")
dim(mirna)
rownames(mirna) <- mirna[,1]
mirna <- mirna[,-1]
exprSet <- mirna[,c(1:4)]
colnames(exprSet) = c("KO1", "KO2", "WT1", "WT2")
exprSet = exprSet[,c(3,4,1,2)]
save(exprSet, file = "./out/file/01_wy_miRNA_counts_matrix.rda")
write.csv(exprSet, file = "./out/file/01_wy_miRNA_counts_matrix.csv")
dat <- exprSet

# 1.Data preprocessing
dat <- t(dat)
dat <- as.data.frame(dat)
dat <- cbind(dat,group_list)
dat[1:4,1:4]
# dat$group_list <- group_list

# 2.Plot PCA diagram
ncol(dat)
dat[,ncol(dat)]
# Drawing
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
plot(dat.pca,choix="ind")

fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = dat$group_list, 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE,
             legend.title = "Groups"
) 

# Variance analysis of DESeq2
class(exprSet)
group <- factor(rep(c("WT", "KO"), times = c(2, 2)), levels = c("WT", "KO"))
group
metadata <- data.frame(sample = c("WT1", "WT2", "KO1", "KO2"), group)
metadata
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = metadata,
                              design = ~group)
nrow(dds)
class(dds)
### Screening genes
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)
### Standardize data using VST methods
### https://mp.weixin.qq.com/s/7BBJGTlOa5i6YPMlrRqw2w
vsd <- varianceStabilizingTransformation(dds, blind = F)

### The assay function extracts the normalized data of VST and saves the data for heat map
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:4,1:4]
save(exprSet_vst, file = "./out/file/03_wy_miRNA_vst.rda")

#################################################################
### DESeq2
dds <- DESeq(dds)
### Use the counts function to extract normalized data or count data
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
save(normalized_counts, file = "./out/file/normalized_counts.Rdata")
head(normalized_counts)

### results Get the results of the variance analysis(KO vs WT)
contrast <- c("group", "KO", "WT")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
### MA
plotMA(dd1, ylim=c(-5,5))
graph2tif(file = "./out/plot/plotMA.tiff")
### logFC correction 
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1, type="normal")
plotMA(dd2, ylim=c(-5,5))
graph2tif(file = "./out/plot/modif_plotMA.tiff")
summary(dd2, alpha = 0.05)
library(dplyr)
library(tibble)
library(tidyr)

### The results of the variance analysis
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("mirna_id") 
class(res)

### Modify the column name
colnames(res) <- c("mirna_id","baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
### Save data
save(res,file = "./out/file/04_wy_mirna_WT_KO_res.rda")
head(res)

# Determine the multiplier of difference
logFC_cutoff <- with(res, mean(abs(logFC)) + 2*sd(abs(logFC)))
logFC_cutoff <- round(logFC_cutoff, 2)
logFC_cutoff
### Screen for differential genes
### Identify up- and down-regulated expression genes
res$change <- as.factor(ifelse(res$P.Value <= 0.05 & abs(res$logFC) > 0.5,
                               ifelse(res$logFC > 0.5, "UP", "DOWN"), "STABLE"))
table(res$change)
write.csv(res, file = "./out/file/wy_mirna_DEG_ALL.csv", row.names = F)
save(res, file = "./out/file/05_wy_mirna_DEG_ALL.rda")




















