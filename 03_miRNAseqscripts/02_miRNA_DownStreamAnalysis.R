### Differential miRNA enrichment analysis
## upmiRNA
library(XML)
library(methods)
anno <- read.table(file = "./data/all_target_gene.xls.annotate.txt", header = T, sep = "\t")
gdes <- read.delim("./resouce/gene.description.xls")
target <- read.delim("./data/all_target_gene.xls")
dndf = res[which(res$change == "DOWN"), ]
dngene = dndf$mirna_id
dngene
length(dngene)
dntarget <- merge(target, dndf, by = "mirna_id")
colnames(dntarget)
dng <- dntarget$target_gene
dns <- list()
for (i in 1:length(dng)) {
  print(i)
  dns[[i]] <- gdes[which(gdes$target_gene == dng[i]), ]
}
dns[[1]]
dnsdf <- do.call(rbind, dns)
head(dnsdf)
dngsdf <- cbind(dntarget, dnsdf[,-1])
dnsymbol <- dngsdf$Gene.name
head(dngsdf)
write.csv(dngsdf, file = "./out/file/09_dn_target_match.csv", row.names = F)

##### Enrichment analysis
library(clusterProfiler)
Ensmus <- unique(dnsymbol)
length(Ensmus)
## GO
#BiocManager::install("org.Mm.eg.db", force = T)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(stringr)
keytypes(org.Mm.eg.db)
gs <- c("Ncf2", "Cybb")
gsf <- bitr(gs, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Mm.eg.db)
gsf

dnentrez <- mapIds(org.Mm.eg.db,
                   keys=Ensmus,
                   column="ENTREZID",
                   keytype="SYMBOL",
                   multiVals="first")

ego_CC <- enrichGO(gene = dnentrez,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2,
                   readable = T)
package.version("clusterProfiler")
dotplot(ego_CC)
write.csv(as.data.frame(ego_CC@result), file = "./out/file/08_wy_mirna_dnego_CC_ALL.csv")
library(export)
graph2tif(file = "./out/plot/08_wy_mirna_dnego_CC.tiff")
graph2pdf(file = "./out/plot/08_wy_mirna_dnego_CC.pdf")

ego_BP <- enrichGO(gene = dnentrez,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2,
                   readable = T)
bpdf <- as.data.frame(ego_BP@result)
write.csv(bpdf, file = "./out/file/08_wy_mirna_dnego_BP_ALL.csv")
dotplot(ego_BP) + scale_y_discrete(labels=function(x) str_wrap(x, width=100))
graph2tif(file = "./out/plot/08_wy_mirna_dnego_BP.tiff")
graph2pdf(file = "./out/plot/08_wy_mirna_dnego_BP.pdf")

ego_MF <- enrichGO(gene = dnentrez,
                   OrgDb = org.Mm.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2,
                   readable = T)
mfdf <- as.data.frame(ego_MF@result)
write.csv(mfdf, file = "./out/file/mirna_dnego_MF_ALL.csv")
dotplot(ego_MF) + scale_y_discrete(labels=function(x) str_wrap(x, width=100))
graph2tif(file = "./out/plot/08_wy_mirna_dnego_MF.tiff")
graph2pdf(file = "./out/plot/08_wy_mirna_dnego_MF.pdf")

## KEGG
# remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
library(KEGG.db)
EGG.dn <- enrichKEGG(gene = dnentrez,
                     organism = 'mmu',
                     pvalueCutoff = 0.2,
                     use_internal_data = T)
KEGG_ddf <- as.data.frame(EGG.dn@result)
eid <- KEGG_ddf["geneID"]
id <- c()
symbol <- list()
genenames <- list()
dim(eid)
options(connectionObserver = NULL)
library("org.Mm.eg.db")
for (i in 1:length(rownames(eid))) {
  print(i)
  id=unlist(strsplit(eid$geneID[i], split = "/"))
  symbol[[i]]=bitr(id, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
  genenames[[i]]=paste(symbol[[i]]$SYMBOL, collapse = "/")
}
gid <- do.call(rbind, genenames)
dim(gid)
eggdf <- cbind(KEGG_ddf, gid)
write.csv(eggdf, file = "./out/file/09_mirna_KEGG_dntargetdf.csv")
browseKEGG(KEGG_udf, "mmu04666")
