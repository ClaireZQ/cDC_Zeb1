### Identifying differential binding sites
remove.packages('rlang')
install.packages("munsell")
install.packages("rlang")
library(tidyverse)
install.packages("parsnip")
BiocManager::install("DiffBind", force = TRUE)

library(DiffBind)
rm(list = ls())
# args <- commandArgs(TRUE)                                                                                   
# samplelist <- "diffbind_samplelist"
out_dir <- "./"
fold <- 0.3
samples <- read.csv("./data/diffbind/Samle_list.txt",sep = "\t")
obj <- dba(sampleSheet=samples)
count <- dba.count(obj, summits=250)
norm <- dba.normalize(count, method=DBA_DESEQ2, library=DBA_LIBSIZE_FULL, normalize=DBA_NORM_LIB)
contrast <- dba.contrast(norm, contrast=c("Condition", "B", "A")) # treat vs control
analyze <- dba.analyze(contrast, method=DBA_DESEQ2)
res <- dba.report(analyze, method=DBA_DESEQ2, contrast = 1, th=1)
save(res, file = "./diffbind_res.rda")
# table
gain <- res[(res$FDR < 0.05) & (res$Fold > fold) ,]
loss <- res[(res$FDR < 0.05) & (res$Fold < -fold) ,]
gain_out <- as.data.frame(gain)
loss_out <- as.data.frame(loss)
f_gain <- sprintf("%s/diffbind_gain.txt", out_dir)
f_loss <- sprintf("%s/diffbind_loss.txt", out_dir)
write.table(gain_out, file=f_gain, sep="\t", quote=FALSE, col.names=NA)
write.table(loss_out, file=f_loss, sep="\t", quote=FALSE, col.names=NA)
# figure
pdf("heat.pdf")
dba.plotHeatmap(norm)
dev.off()


class(res)
length(res@seqnames)
head(res@seqnames)
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(GenomicFeatures))
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(clusterProfiler)
# rm(list = ls())
# gff <- args[1]
# peakFile <- args[2]
# name <- args[3]

# tx <- makeTxDbFromGFF(file="./resouce/mm10.refGene.gtf")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# peak <- readPeakFile(peakfile="./data/cDC1c_peaks.narrowPeak",header=FALSE,as="GRanges")
d <- mcols(res) # input Diffbind result
if(ncol(d)==1){
  colnames(d)="peakName"
}else if (ncol(d)==7){
  colnames(d)=c("peakName","colorScore","strand","fold_enrichment","-log10(pvalue)","-log10(qvalue)","relative_summit_position")
  d=d[,-3]
}

mcols(res) <- d

peakAnno_diffbind <- annotatePeak(peak=res, tssRegion = c(-3000, 3000), 
                         TxDb=txdb, annoDb = "org.Mm.eg.db",
                         assignGenomicAnnotation=TRUE)

peakAnno.result <- as.data.frame(peakAnno_diffbind@anno)
# table
fold = 0.3
gain <- peakAnno.result[(peakAnno.result$p.value < 0.05) & (peakAnno.result$Fold > fold) ,]
loss <- peakAnno.result[(peakAnno.result$p.value < 0.05) & (peakAnno.result$Fold < -fold) ,]
gain_out <- as.data.frame(gain)
loss_out <- as.data.frame(loss)
f_gain1 <- sprintf("%s/diffbind_gain.txt", out_dir)
f_loss1 <- sprintf("%s/diffbind_loss.txt", out_dir)
write.table(gain_out, file=f_gain, sep="\t", quote=FALSE, col.names=NA)
write.table(loss_out, file=f_loss, sep="\t", quote=FALSE, col.names=NA)




write.table(peakAnno.result,
            file=paste0("diffbind_DEPeak", ".peakAnno.txt"),
            sep="\t",
            row.names=FALSE)
















