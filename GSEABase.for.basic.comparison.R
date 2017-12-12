#GSEABase: GSEA analysis on DE genes per condition
#October 14th, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/GO_analysis_by_DAVID_Oct_2015/")
library(GSEABase)

#Create an expression dataset
list.of.name.of.conditions <- list("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d",
                                   "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                   "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
                                   "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
normalized.counts = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_mega_RNA-seq_Sept_2015/Mega_RNA-seq_edgeR_normalized_counts_with_more_than_5_counts_averaged_with_gene_name_and_ID.txt", header=T)

for (i in 1:length(list.of.name.of.conditions)){
    table = read.table(file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/GO_analysis_by_DAVID_Oct_2015/list_of_DE_genes_in_", list.of.name.of.conditions[i],".txt", sep=""), header=T)
    list.of.genes = table[,1]
    table.with.values = data.frame(ncol=3, nrow=length(list.of.genes))
    for (j in 1:length(list.of.genes)){
        table.with.values[j,1] = as.character(normalized.counts[which(normalized.counts$gene_id == as.character(list.of.genes[j])),1])
        table.with.values[j,2] = as.character(normalized.counts[which(normalized.counts$gene_id == as.character(list.of.genes[j])),2])
        table.with.values[j,3] = as.numeric(normalized.counts[which(normalized.counts$gene_id == as.character(list.of.genes[j])),i+2])
    }
    colnames(table.with.values) = c("gene_id", "gene_name", "nor.expre.value")
    write.table(table.with.values, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/GSEABase_prep/", list.of.name.of.conditions[i],"_DE_genes_with_exp_values.txt", sep=""), quote=F, col.names=T, row.names=F)
}