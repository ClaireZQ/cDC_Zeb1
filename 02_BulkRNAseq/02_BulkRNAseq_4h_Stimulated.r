############################################################################
############################################################################
### Author：QuanZhang
### time：2023-05-05
### Data analysis of stimulated RNAseq for 4h
############################################################################

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

getwd()
dir.create("./new20220219")
rm(list=ls())
## read data
expr <- data.table::fread(file = "./data/gene_count.csv")
class(expr)
expr <- as.data.frame(expr)
class(expr)
dim(expr)
expr[1:3,1:3]
data <- expr[,c(1:5)]
head(data)
colnames(data)
colnames(data)[2:5]=c("KO1","WT2","KO2","WT1")
head(data)
data=data[,c(1,5,3,2,4)]
head(data)
# Check matrix
library(dplyr)
library(tidyr)
library(tibble)
dim(data)
head(data)
exprSet <- data %>% 
  column_to_rownames("gene_id") 
# View low/zero expression genes
table(rowSums(exprSet==0)==ncol(exprSet))
# Remove low- and non-expressed genes
exprpca <- exprSet[apply(exprSet, 1, sum) != 0, ]

# Using the original expression data for the next processing
library(edgeR)
exprpca <- cpm(exprpca, log = T)

# The absolute median difference MAD/standard deviation SD statistical method was used for data outlier detection 
tmp <- sort(apply(exprpca,1, mad),decreasing = T)[1:500]
exprpca <- exprpca[names(tmp),]
# The expression of 500 genes was used to make a correlation map
library(corrplot)
dim(exprpca)
# Calculate relevance
M <- cor(exprpca)
g <- corrplot(M,method = "color")
# Heat map sample correlation
pheatmap::pheatmap(M)

# Scale plots heat maps after standardizing variables
dir.create("./new20220219/plot")
pheatmap::pheatmap(scale(M))
library(export)
graph2tif(file = "./new20220219/plot/01_sample_corrplot.tiff")
graph2pdf(file = "./new20220219/plot/01_sample_corrplot.pdf")

dir.create("./new20220219/file")
save(exprSet, data, file = "./new20220219/file/00_expression.rda")

# Check data
group_list <- c(rep("WT", 2), rep("KO", 2))
group_list
exprSet[1:4,1:4]
dim(exprSet)
dat <- exprSet
# 1.Data preprocessing
dat <- t(dat)
dat <- as.data.frame(dat)
dat <- cbind(dat,group_list)
dat[1:4,1:6]

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
library(DESeq2)
group <- factor(rep(c("WT", "KO"), times = c(2, 2)), levels = c("WT", "KO"))
group
metadata <- data.frame(sample=c("WT1","WT2","KO1","KO2"), group)
class(metadata$sample)
class(metadata$group)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ group,
                              tidy = TRUE)

nrow(dds)
### Screening genes
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)
### Standardize data using VST methods
### https://mp.weixin.qq.com/s/7BBJGTlOa5i6YPMlrRqw2w
vsd <- vst(dds, blind = FALSE)
class(dds)

# Save the preprocessed expression matrix
save(vsd, file = "./new20220219/file/01_wy_4h_exprSet_vsd.rda")

### The assay function extracts the normalized data of VST and saves the data for heat map
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:4,1:4]
save(exprSet_vst, file = "./new20220219/file/02_HKLM_wy_4h_exprSet_vst.rda")

#################################################################
### DESeq2
dds <- DESeq(dds)
### Use the counts function to extract normalized data or count data
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
save(normalized_counts,file = "./new20220219/file/03_HKLM_wy_4h_normalized_counts.rda")
head(normalized_counts)
### results Get the results of the variance analysis(KO vs WT)
contrast <- c("group", "KO", "WT")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
### MA
plotMA(dd1, ylim=c(-5,5))
graph2tif(file = "./new20220219/plot/plot_MA.tiff")
graph2pdf(file = "./new20220219/plot/plot_MA.pdf")
### logFC correction 
dd2 <- lfcShrink(dds, contrast=contrast, res=dd1,type="normal")
plotMA(dd2, ylim=c(-5,5))
graph2tif(file = "./new20220219/plot/plot_MA_ifcshrink.tiff")
graph2pdf(file = "./new20220219/plot/plot_MA_ifcshrink.pdf")
summary(dd2, alpha = 0.05)

### The results of the variance analysis
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") #%>% 
  #separate(gene_id,into = c("gene_id"),sep = "\\.") 
class(res)

#### Annotation
rest <- merge(res, expr, by = "gene_id")
colnames(rest)
rest$entrez <- mapIds(org.Mm.eg.db,
                      keys=rest$gene,
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
colnames(rest)
### Modify the column name
colnames(rest)[c(1:7, 12)] <- c("gene_id","baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val", "gene")
colnames(rest)

### Save data
save(rest, res, file = "./new20220219/file/04_HKLM_wy_4h_DESeq2_rest_res.rda")

# Determine the multiplier of difference
logFC_cutoff <- with(rest, mean(abs(logFC)) + 2*sd(abs(logFC)))
logFC_cutoff <- round(logFC_cutoff, 2)
logFC_cutoff

### Screen for differential genes
### Identify up- and down-regulated expression genes
rm(list = ls())
rest$change <- as.factor(ifelse(rest$P.Value < 0.05 & abs(rest$logFC) > 0.6,
                                ifelse(rest$logFC > 0.6, "UP", "DOWN"), "STABLE"))


table(rest$change)
write.csv(rest, file = "./new20220219/file/05_HKLM_4h_wy_deg_all_Fc3.csv", row.names = F)
save(rest, file = "./new20220219/file/05_HKLM_4h_wy_deg_all_Fc3.rda")

#### FC1.5&Pvalue0.05
rm(list = ls())
rest <- read.csv(file = "./new20220219/file/05_HKLM_4h_wy_deg_all.csv",header = T, row.names = 1,sep = ",")
load(file = "./new20220219/file/05_HKLM_4h_wy_deg_all.rda")
table(rest$change)


#########################################################################
### GO 分析
library(clusterProfiler)
library(stringr)
deg=rest[which(rest$change=="DOWN"|rest$change=="UP"),]
dim(deg)
save(deg, file = "./new20220219/file/06_HMKL_wy_4h_deg_up&down.rda")

gene_df <- deg %>% 
  dplyr::select(gene,logFC,entrez, change) %>% 
  filter(entrez!="NA") %>% 
  distinct(entrez,.keep_all = T)
  
ego <- enrichGO(gene = gene_df$entrez,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                minGSSize = 1,
                readable = TRUE)
save(ego, file = "./new20220219/file/07_HKLM_wy_4h_ego_degALL.rda")
goo <- data.frame(Category = ego$ONTOLOGY,
                 ID = ego$ID,
                 Term = ego$Description,
                 Genes = gsub("/", ", ", ego$geneID),
                 p.val = ego$pvalue,
                 adj_pval = ego$p.adjust)
class(goo)
write.csv(goo, file = "./new20220219/file/07_HKLM_wy_4h_ego_degALL.csv")


### KEGG
library(KEGG.db)
kegg <- enrichKEGG(gene = gene_df$entrez,
                   keyType = "kegg",
                   organism = "mmu",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.2,
                   pAdjustMethod = "BH",
                   use_internal_data = F)
barplot(kegg)
dotplot(kegg)
head(kegg@result)
save(kegg, file = "./new20220219/file/08_HKLM_wy_4h_kegg.rda")

#########################################################################
### GSEA
rm(list = ls())
library(msigdbr)
library(clusterProfiler)
packageVersion("msigdbr")
msigdbr_species()
library(dplyr)
mm=msigdbr(species = "Mus musculus")
mm %>%
  distinct(gs_cat, gs_subcat) %>%
  arrange(gs_cat, gs_subcat)

gene_df <- rest %>%
  dplyr::select(gene_id, logFC, gene, entrez) %>%
  filter(entrez != "NA") %>%
  distinct(entrez, .keep_all = T)
genelist <- gene_df$logFC
names(genelist) = gene_df$gene
genelist = sort(genelist, decreasing = T)
head(genelist)
m_df=msigdbr(species = "Mus musculus", category = "H")
mhallmarks <- m_df %>%
  dplyr::select(gs_name,gene_symbol)
m_dfk=msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
m_dfR=msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
class(mhallmarks)
mkegg <- m_dfk %>%
  dplyr::select(gs_name, gene_symbol)
mreac <- m_dfR %>%
  dplyr::select(gs_name, gene_symbol)
mhall <- m_df %>%
  dplyr::select(gs_name, gene_symbol)

allgsea <- GSEA(genelist, TERM2GENE = mkegg)
allgd=as.data.frame(allgsea@result)
write.csv(allgd, file = "./new20220219/file/16_kegg_gsea.csv")










