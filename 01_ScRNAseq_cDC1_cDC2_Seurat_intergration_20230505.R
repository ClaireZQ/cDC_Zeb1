####load seuratb packages
library(Seurat)
library(dplyr)
library(umap)
library(ggplot2)
library(cowplot)
library(xlsx)
load('/home/zhangquan/project/singlecell/BDxiao/DC-01/2000features/Seurat_object_DC1.RData')
load('/home/zhangquan/project/singlecell/BDxiao/DC-02/2000features/Seurat_object_DC2.RData')
####read DC1 data and normalize
DC1
DC1@meta.data$sample <- 'DCWT'
DC1[['batch']] <- DC1[['sample']]
DC1[['percent.mt']] <- PercentageFeatureSet(DC1,pattern = '^mt-')
DC1 <- subset(x = DC1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
DC1 <- NormalizeData(object = DC1, normalization.method = 'LogNormalize',scale.factor = 10000)
DC1 <- FindVariableFeatures(DC1,selection.method = 'vst',nfeatures = 2000)
top10_3 <- head(VariableFeatures(DC1),10)
DC1 <- CellCycleScoring(DC1,s.features = s.genes,g2m.features = g2m.genes,set.ident = T)


####read DC2 data and normalize
DC2
DC2@meta.data$sample <- 'DCKO'
DC2[['batch']] <- DC2[['sample']]
DC2[['percent.mt']] <- PercentageFeatureSet(DC2,pattern = '^mt-')
DC2 <- subset(x = DC2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
DC2 <- NormalizeData(object = DC2, normalization.method = 'LogNormalize',scale.factor = 10000)
DC2 <- FindVariableFeatures(DC2, selection.method = 'vst',nfeatures = 2000)
top10_4 <- head(VariableFeatures(DC2),10)
DC2 <- CellCycleScoring(DC2,s.features = s.genes,g2m.features = g2m.genes,set.ident = T)

####Mulitple setdatas deposited seuratlist
DC.list <- list(DC1,DC2)

####NormalizeData and find different genes
#for (i in 1:length(pancreas.list)) {

#}

####Cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for (i in 1:length(DC.list)) {
DC.list[[i]] <- CellCycleScoring(DC.list[[i]],s.features = s.genes,g2m.features = g2m.genes,set.ident = T)
}

#for (i in 1:length(pancreas.list)) {
#pancreas.list[[i]] <- ScaleData(pancreas.list[[i]],vars.to.regress = c('nCount_RNA','percent.mt','S.Score','G2M.Score'),features = rownames(pancreas.list[[i]]))
#}
#### SCTransform
####for (i in 1:length(pancreas.list)) {
####pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], vars.to.regress = c('nCount_RNA','percent.mt','S.Score','G2M.Score'),verbose = F )
####}

#### PrepSCTIntegration
#pancreas.features <- selectIntegrationFeatures(object.list = pancreas.list,nfeatures = 3000)
#merged@meta.data$clusters <- merged@meta.data$seurat_clusters
#length(merged@meta.data$sample[merged@meta.data$sample == 'DC2'])#3608
#3length(merged@meta.data$sample[merged@meta.data$sample == 'DC1' & merged@meta.data$clusters == 8])#335
#counts1 <- list()
#counts2 <- list()
#for (i in 0:15) {
#	counts1[[i+1]] <- length(rownames(merged@meta.data)[merged@meta.data$seurat_clusters == i & merged@meta.data$sample == 'DC1'])
#	counts2[[i+1]] <- length(rownames(merged@meta.data)[merged@meta.data$seurat_clusters == i & merged@meta.data$sample == 'DC2'])
 
#}

####Find anchor by findintergrationAnchors
DC.anchors <- FindIntegrationAnchors(object.list = Endo, anchor.features = 2000,
dims = 1:30)

####Use these anchors to correct technical and batch differences
DC.integrated <- IntegrateData(anchorset = DC.anchors, dims = 1:30)
merged <- DC.integrated
####Two expression matrices can be switched through DefaultAssay function
DefaultAssay(merged) <- 'integrated'

####Seutat analysis procedure
merged <- ScaleData(merged, vars.to.regress = c('nCount_RNA','percent.mt','S.Score','G2M.Score','batch'),features = rownames(merged))
merged <- RunPCA(merged,npcs=50,verbose = F)
merged <- JackStraw(object = merged,reduction = 'pca',dims = 50,num.replicate = 100)
merged <- ScoreJackStraw(object = merged,dims = 1:50,reduction = 'pca',do.plot = F)
pdf('jackstrawplot.pdf',width = 15,height = 8)
JackStrawPlot(object = merged,dims = 1:50)
dev.off()
pdf('ElbowPlot.pdf',width=8,height=8)
ElbowPlot(object = merged,ndims = 50,reduction = 'pca')
dev.off()


####Dimensional Reduction(UMAP) and Visualization
merged <- RunUMAP(merged,reduction = 'pca', dims = 1:20)
merged <- FindNeighbors(object = merged, reduction = 'pca',dims = 1:20)
merged <- FindClusters(object = merged,resolution = 0.5)
pdf('integrated_Endo1.pdf',width = 15, height = 8)
p1 <- DimPlot(merged, reduction = 'umap', group.by = 'location')
p2 <- DimPlot(merged, reduction = 'umap', label = T)
plot_grid(p1,p2)
dev.off()


####To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
pdf('Endo_splitgraph.pdf',width=15,height=8)
DimPlot(merged,reduction = 'umap', split.by = 'sample',label = T, label.size = 3, pt.size =0.1)
dev.off()

####Identify conserved cell type markers
DefaultAssay(merged) <- 'RNA'
con.markers <- list()
con.markers1 <- list()
con0.markers <- FindConservedMarkers(merged, ident.1 = 0, grouping.var = 'sample', verbose = F)
con8.markers <- FindConservedMarkers(merged, ident.1 = 8, grouping.var = 'sample', verbose = F)
con.markers12 <- FindConservedMarkers(merged, ident.1 = 12, grouping.var = 'sample', verbose = F)
for (i in 18:20) {
con.markers[[i+1]] <- FindConservedMarkers(merged, ident.1 = i, grouping.var = 'sample', verbose = F)
}
names(con.markers) <- c('con18', 'con19')
names(con.markers) <- c('con0','con1','con2','con3','con4','con5','con6','con7','con8','con9','con10')
for (i in 1:11) {write.xlsx(con.markers[[i]],file='conmarker0_10.xlsx',sheetName=names(con.markers)[i],
                     col.names=T,row.names=T,append=T)}
for (i in 1:3) {write.xlsx(con.markers[[i]],file='conmarker1_3.xlsx',sheetName=names(con.markers)[i],
                     col.names=T,row.names=T,append=T)}
for (i in 19:20) {write.xlsx(con.markers[[i]],file='conmarker20_19.xlsx',sheetName=names(con.markers)[i],
                     col.names=T,row.names=T,append=T)}
for (i in 14:15) {write.xlsx(con.markers[[i]],file='conmarker13_14.xlsx',sheetName=names(con.markers)[i],
                     col.names=T,row.names=T,append=T)}


write.xlsx(con.markers[[20]],file='conmarker19.xlsx',sheetName='con19',
                     col.names=T,row.names=T,append=T)



#we can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
pdf('annotate_features_clusters3_4.pdf',width = 10, height = 10)
FeaturePlot(merged, features =  c('Cd24a','Siglech','Ly6c2','Ccr2','Ly6d','Ifi205','Anxa1','Id2','Aif1'),min.cutoff = 'q9')
dev.off()

pdf('VDR.pdf',width = 10, height = 10)
FeaturePlot(merged, features = 'VDR', min.cutoff = 'q1')
dev.off()


####renames clusters
levels(x=merged)
#DC1_2
merged <- RenameIdents(merged, '0' = 'cDC2', '1' = 'cDC2', '2' = 'cDC2', '3' = 'cDC2', '4' = 'cDC1',
	'5' = 'cDC2', '6' = 'cDC2', '7' = 'cDC1', '8' = 'cDC2', '9' = 'cDC2', '10' = 'cDC1', '11' = 'cDC2',
	'12' = 'cDC2', '13' = 'cDC2', '14' = 'cDC2', '15' = 'cDC2', '16' = 'cDC2', '17' = 'cDC2')
merged <- RenameIdents(merged, 'cDC2' = '0', 'cDC2' = '1', 'cDC2' = '2', 'cDC2' = '3', 'cDC1' = '4',
	'cDC2' = '5','cDC2' = '6', 'cDC1'='7', 'cDC2' = '8','cDC2'='9', 'cDC1'='10', 'cDC2'='11',
	'cDC2'='12', 'cDC2'='13', 'cDC2'='14', 'cDC2'='15', 'cDC2'='16', 'cDC2'='17')
head(x = Idents(object = merged))

merged@meta.data$cell_ann = merged@meta.data$seurat_clusters
merged@meta.data$cell_ann[merged@meta.data$cell_ann=='4'|merged@meta.data$cell_ann=='7']='cDC1'

merged@meta.data$type[1:2581]='WT'
merged@meta.data$type[2582:7101]='KO'
#if(merged@meta.data$cell_ann != 'cDC1'){merged@meta.data$cell_ann=='cDC2'}else{merged@meta.data$cell_ann=='cDC1'}
merged@meta.data[which(merged@meta.data$cell_ann != 'cDC1'),'cell_ann']='cDC2'
merged@meta.data$cell_type=paste(merged@meta.data$cell_ann,merged@meta.data$type,sep = '_')
pdf('annotate_cellclusters_DC1_2_split_ident.pdf',width = 10,height = 8)
DimPlot(merged, label = T, size = 8, split.by = 'ident')
dev.off()
pdf('annotate_cellclusters_DC1_2_split.pdf',width = 18,height = 10)
DimPlot(merged, reduction = 'umap',label = T, size = 8, split.by = 'sample',group.by = 'cell_ann')
dev.off()

pdf('annotate_cellclusters_DC1_2_merged.pdf',width = 15,height = 12)
DimPlot(merged, label = T, size = 8)
dev.off()

#### Comparison of genes in two different sets of clusters
#merged@meta.data$celltype <- c(rep('WT',3608), rep('KO',6032))
merged@meta.data$celltype <- paste(Idents(merged),merged@meta.data$sample,sep = '_')
merged$cluster <- Idents(merged)
Idents(merged) <- 'celltype'

####split.by genes
pdf('genes_differen_clusters3_4.pdf',width = 8, height = 40)
FeaturePlot(merged,features = c('Cd24a','Siglech','Ly6c2','Ccr2','Ly6d','Ifi205','Anxa1','Id2','Aif1'),
split.by = 'sample', max.cutoff = 9, cols = c('grey','red'))
dev.off()

pdf('genes_xcr1.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Xcr1',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Sirpa.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Sirpa',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Mapk13.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Mapk13',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Tnfsf9.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Tnfsf9',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Tnf.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Tnf',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Irf4.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Irf4',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_Irf8.pdf',width =8, height=4)
FeaturePlot(merged,features = 'Irf8',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()

pdf('genes_B2m.pdf',width =8, height=4)
FeaturePlot(merged,features = 'B2m',
split.by = 'sample', max.cutoff = 1, cols = c('grey','red'))
dev.off()



####VlnPlot
pdf('vln_merged_genes_split3_4.pdf',width = 12, height = 50)
plots <- VlnPlot(merged, features = c('Cd24a','Siglech','Ly6c2','Ccr2','Ly6d','Ifi205','Anxa1','Id2','Aif1'),split.by = 'sample',
group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_apol7c.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Apol7c',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

library(patchwork)
pdf('genesvln_Xcr1.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Xcr1',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_Sirpa.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Sirpa',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_IRF4.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Irf4',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_IRF8.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Irf8',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_B2m.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'B2m',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf('genesvln_tnf.pdf',width =8, height=4)
plots <- VlnPlot(merged,features = 'Tnf',
split.by = 'sample',group.by = 'cluster', pt.size = 0, combine = F)
wrap_plots(plots = plots, ncol = 1)
dev.off()

####Extract original cell expression of the genes 
integrated.counts <- merged@assays$integrated['Xcr1',]
integrated.counts <- as.data.frame(integrated.counts)
sum(integrated.counts == 0)


####Fisher test to calculate the significance proportion of cells with positive expressed genes in the two groups
gene1 <- DC1@assays$RNA[c('Xcr1','Zeb1','Zeb2','Irf8','Ahr','Itgae','Itgam'),]
gene1 <- as.data.frame(gene1)
gene2 <- DC2@assays$RNA[c('Xcr1','Zeb1','Zeb2','Irf8','Ahr','Itgae','Itgam','Clec9a','Batf3','Klf4','Notch2','Irf4','Sirpa','Esam'),]
gene2 <- as.data.frame(gene2)
cell1 <- ncol(gene1)
cell2 <- ncol(gene2)
fish_result <- NULL
for (i in 1:nrow(gene1)) {
print(i)
a <- sum(gene1[i,]==0)
b <- sum(gene2[i,]==0)
df <- matrix(c(cell1-a,cell2-b,a,b), 2, 2, byrow = T)
fish_result1 <- fisher.test(df)	
pvalue=fish_result1$p.value
OR=fish_result1$estimate
CIlow=fish_result1$conf.int[1]
CIhigh=fish_result1$conf.int[2]
genes=rownames(gene1)[i]
fish_result1 <- data.frame(genes,pvalue, OR, CIlow, CIhigh)
fish_result <- rbind(fish_result,fish_result1)
}
write.table(fish_result,file = '14genesfishtest3.txt',row.names=F,sep='\t')

####Use ggplot to draw forestplot
pdf('genesforestp1.pdf',width = 8, height = 6)
p <- ggplot(fish_result, aes(x=genes, y=OR, ymin = CIlow, ymax= CIhigh, colour=genes))+
geom_pointrange(shape=15,size = 0.8,position=position_dodge(width=c(0.1)))+
coord_flip() +
xlab('genes') +
theme_bw() 
for(i in 1:nrow(fish_result)){
	p <- p+annotate('text',x = fish_result$genes[i], y = fish_result$OR[i]+0.3, label=signif(fish_result$pvalue[i],3))
}
print(p)
dev.off()

####Draw scatter plot
wtxcr1 <- merged@assays$RNA['Xcr1',]
wtxcr1 <- as.data.frame(wtxcr1)
wtxcr1 <- wtxcr1[,1:3608]
wtintexcr1 <- merged@assays$integrated['Xcr1',]
wtintexcr1 <- as.data.frame(wtintexcr1)
wtintexcr1 <- wtintexcr1[,1:3608]
row.names(wtintexcr1) <- 'con_Xcr1'
wtxcr1$gene <- 'Xcr1'
wtintexcr1$gene <- 'con_Xcr1'
dfxcr1 <- rbind(wtxcr1,wtintexcr1)
dfxcr1_2 <- t(dfxcr1_2)
dfxcr1_2 <- as.data.frame(dfxcr1_2)
df <- melt(dfxcr1_2)
pdf('Xcr1_disbution2.pdf')
ggboxplot(df, x = 'variable', y = 'value', color = 'variable',
palette = c('#00AFBB','#FC4E07'),
add = 'jitter', shape = 'variable' )
dev.off()

pdf('Xcr1_disbution1.pdf')
ggviolin(df, x = 'variable', y = 'value', fill = 'variable',
palette = c('#00AFBB','#FC4E07'),
add = 'boxplot', add.params = list(fill = 'white') )
dev.off()

expr <- df$value
my_df <- data.frame(Xcr1=expr[1:3608], con_Xcr1=expr[3609:9640])
pdf('scatter.pdf')
ggplot(my_df, aes(x = Xcr1, y = con_Xcr1))+
geom_point()
dev.off()

pdf('scatter2.pdf')
ggplot(my_df, aes(x = Xcr1, y = con_Xcr1)) +
geom_point(shape = 19)+
xlab('Xcr1') + ylab('con_Xcr1') +
geom_smooth(method = lm)
dev.off()

####Identify differential expressed genes across conditions
theme_set(theme_cowplot())
Idents(merged) <- 'celltype'
cl4.cells <- subset(merged,idents = 'cluster4')
Idents(cl4.cells) <- 'sample'
avg.cl4.cells <- log1p(AverageExpression(cl4.cells,verbose=F)$RNA)
avg.cl4.cells$gene <- rownames(avg.cl4.cells)

cl7.cells <- subset(merged,idents = 'cluster7')
Idents(cl7.cells) <- 'sample'
avg.cl7.cells <- log1p(AverageExpression(cl7.cells,verbose=F)$RNA)
avg.cl7.cells$gene <- rownames(avg.cl7.cells)


genes.to.label = c('Apol7c','Xcr1','Zeb1','Zeb2','Irf8','Ahr','Itgae','Itgam','Clec9a','Batf3','Id2','Klf4','Notch2','Irf4','Sirpa','Esam')
pdf('cl4_cl7_identify_deg_acrossconditions.pdf',width = 8, height=8)
p1 <- ggplot(avg.cl4.cells,aes(DCWT,DCKO))+
geom_point()+
ggtitle('cluster4.cells')
p1 <- LabelPoints(plot=p1, points = genes.to.label, repel = T)
p2 <- ggplot(avg.cl7.cells, aes(DCWT,DCKO))+
geom_point()+
ggtitle('cluster7.cells')
p2 <- LabelPoints(plot=p2, points = genes.to.label, repel = T)
plot_grid(p1,p2)
dev.off()

#Kegg and GO analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(xlsx)
library(stringr)
library(readxl)
keytypes(org.Mm.eg.db)

sheet.index <- c(1:18)
data.list <- list()
for (i in sheet.index) {
	data.list[[i]] <- read_excel('kwdemarker.xlsx',sheet=i,col_names=T)
}
for (i in sheet.index) {
	colnames(data.list[[i]])[1] = 'gene'
}
for (i in sheet.index) {
	data.list[[i]] <- as.data.frame(data.list[[i]])
}
gene.df <- list() 
for (i in sheet.index) {
	gene.df[[i]] <- bitr(data.list[[i]][,'gene'], fromType = "SYMBOL", 
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db)
}

ggo <- lapply(c(1:18), function(x){
	enrichGO(gene           = gene.df[[x]]$ENTREZID,
		     OrgDb          = org.Mm.eg.db,
		     ont            = 'ALL',
		     pAdjustMethod  = 'BH',
		     qvalueCutoff   = 0.05,
		     readable       = T)
})

for (i in 16:18) {
	write.xlsx(ggo[[i]],file = 'ggo1.xlsx',sheetName=names(ggo)[i], col.names = T, row.names = T, append = T)
}	
write.xlsx(ggo[[15]],file = 'ggo1.xlsx',sheetName=names(ggo)[15], col.names = T, row.names = T, append = T)
write.csv(ggo[[15]],file = 'ggo15.csv',quote=F)
write.csv(ggo[[16]],file = 'ggo16.csv',quote=F)
write.csv(ggo[[17]],file = 'ggo17.csv',quote=F)
write.csv(ggo[[18]],file = 'ggo18.csv',quote=F)


kegg.result <- lapply(c(1:18), function(x){
       enrichKEGG(unique(gene.df[[x]][,3]),
             	  organism = 'mmu', 
	              keyType = 'kegg', 
	              pvalueCutoff = 0.05,
	              pAdjustMethod = 'fdr', 
	              qvalueCutoff = 0.1,
	              minGSSize = 10,
	              maxGSSize = 300,
	              use_internal_data = FALSE)
})


for(i in 1:18){
	kegg.result[[i]]@result$geneName <- NULL
	#kegg.result[[i]]@result$TF <- NULL
	for(j in 1:nrow(kegg.result[[i]]@result)){

		s <- str_split_fixed(kegg.result[[i]]@result[j,'geneID'], pattern='/', n=100)
    	s <- s[1:length(which(s!=""))]

    	genetmp <- bitr(s, fromType = "ENTREZID", 
       				toType = c("SYMBOL"),
       				OrgDb = org.Mm.eg.db)
    	kegg.result[[i]]@result[j,'geneName'] <- paste(genetmp[,'SYMBOL'],collapse = ", ")

	}
	
}

cluster <- paste0('cluster',0:17)		
names(kegg.result) <- cluster
for(i in 1:18){ write.xlsx(kegg.result[[i]]@result,file="Allgene_kegg_all.xlsx",sheetName=names(kegg.result)[i],col.names=T,row.names=T,append=T)}



