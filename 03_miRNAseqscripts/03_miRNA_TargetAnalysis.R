##### Use different web to predict target genes
#### Target genes are predicted by ENCORI
### First download the table of the regulatory relationship between miRNA and mRNA
library(stringr)
# family <- read.table("./resouce/starBase_refData/mm10_all_fimaly.txt", sep = "\t")
# class(family)
# miRNA <- unlist(strsplit(as.character(family$V4), ","))
# length(miRNA)
# dir.create("./out/starBasemRNA")
# for (mir in miRNA) {
#   file = paste("./out/starBasemRNA/", mir, ".txt", sep = "")
#   link=paste("http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=mm10&geneType=mRNA&miRNA=",mir,"&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all", sep = "")
#   download.file(link, file)
#   Sys.sleep(5)
# }
###########################################
#combine mRNA_miRNA_interaction
###########################################
# mRNA_files <- list.files(path="mRNA", full.names=TRUE)
# 
# library(plyr)
# lncRNA.list <- llply(mRNA_files, function(x)read.table(x,header=T,sep="\t",comment.char ="#",stringsAsFactors=F))
# combind_mRNA=do.call(rbind,mRNA.list)
# write.table(file="mRNA_miRNA_interaction.txt",combind_mRNA,quote=F,sep="\t",row.names = F)

##########################################
# miRNA-mRNA，mRNA-miRNA
##########################################
# Read the regulatory relationship data between the original mRNA and the miRNA
mir_mRNA=read.table("./resouce/mRNA_miRNA_interaction.txt", sep="\t", stringsAsFactors=F, header=T)
mir_mRNA=unique(mir_mRNA[,c("miRNAname","geneName")])
head(mir_mRNA)
# Convert the data frame to list, with the name of the gene named list
mRNA_mir_list=unstack(mir_mRNA,miRNAname~geneName)
# Convert the data frame to a list, with the miRNA name as the name of the list
mir_mRNA_list=unstack(mir_mRNA,geneName~miRNAname)

### upgene
mir_candidate=upgene[1:16]
# Get a subset of the list, you can get the regulatory relationship between miRNA-mRNA, and then convert it into a data frame
mir_target_mRNA=stack(mir_mRNA_list[mir_candidate])
# Change the columnnames of the data frame
names(mir_target_mRNA)=c("gene","miRNA")
dim(mir_target_mRNA)
table(mir_target_mRNA)
length(unique(mir_target_mRNA$miRNA))
smirna <- as.character(unique(mir_target_mRNA$miRNA))
mirna_target_df <- list()
convert <- function(smirna){
  gid <- mir_target_mRNA[which(mir_target_mRNA$miRNA == smirna), "gene"]
  symbol = paste(gid[1:length(gid)], collapse = ",")
  mirna_target_df <- data.frame(miRNA = smirna, symbol)
  }
mi_ta_df <- do.call(rbind, lapply(smirna, convert))
# Write out the knot to a file
write.table(file="./out/file/starbase_miRNA_to_mRNA.txt",mir_target_mRNA,quote=F,sep="\t",row.names = F)
write.table(file = "./out/file/starbase_miRNA_to_mRNA_convert.xls", mi_ta_df, quote = F, sep = "\t", row.names = F)



#### Target genes are predicted by TargetScan
############################################
# mir to gene 
###########################################
# Read miRNA family information
mir_family=read.table("./resouce/targetscan/miR_Family_Info.txt",header=T,sep="\t",stringsAsFactors = F)
# The miRNA family information of humans was extracted, and the species ID number of humans was 9606 and that of rats was 10090
mouse_mir_family=mir_family[mir_family$Species.ID==10090,c("miR.family","MiRBase.ID")]

# Read all predictions for Targetscan
predicted_targets=read.table("./resouce/targetscan/Predicted_Targets_Info.default_predictions.txt",header=T,sep="\t",stringsAsFactors = F)
# Extract human predictions based on species number
mouse_targets=predicted_targets[predicted_targets$Species.ID==10090,]
# Convert the regulatory relationship data.frame into a list, whose name is the name of the miRNA family
mouse_targets_list=unstack(mouse_targets,Gene.Symbol~miR.Family)

 
mir_candidate=upmirna[1:16]
# Obtain the miRNA family corresponding to the miRNA
index=match(mir_candidate,mouse_mir_family$MiRBase.ID)
mir_candidate_family=mouse_mir_family[index,"miR.family"]

# A subset of the list was extracted to obtain the target gene of the miRNA family
mir_targets=mouse_targets_list[mir_candidate_family]
# Replace the name of the list with the name of miRNA
names(mir_targets)=mir_candidate


### list to dataframe 
mir_targets_table=stack(mir_targets)
names(mir_targets_table)=c("gene","mir")

length(unique(mir_targets_table$mir))
tmirna <- as.character(unique(mir_targets_table$mir))
tmirna_target_df <- list()
tconvert <- function(tmirna){
  gid <- mir_targets_table[which(mir_targets_table$mir == tmirna), "gene"]
  symbol = paste(gid[1:length(gid)], collapse = ",")
  tmirna_target_df <- data.frame(miRNA = tmirna, symbol)
}
tmi_ta_df <- do.call(rbind, lapply(tmirna, tconvert))
# Write out the knot to a file
write.table(file="./out/file/targetscan_mir_to_targets.txt",mir_targets_table,quote=F,sep="\t",row.names=F)
write.table(file="./out/file/targetscan_mir_to_targets_convert.xls",tmi_ta_df,quote=F,sep="\t",row.names=F)
