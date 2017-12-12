#Answering Nicolas' questions
#April 18th 2017
#Joo Hyun Im (ji72)

####################################################
#1. Identify the genes in each PC
####################################################

rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq")

#From /Users/JooHyun/Dropbox/Cornell/Lab/R.scripts/deseq.with.time.R (Dec 2014/Jan 2015)
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)

#for 12hr 10 live infections + heatkilled
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103),]
#for 12hr 10 live infections only
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99),]
#for 12hr 9 live infections except Staph + heatkilled
CountTable <- CountTable[,c(8,11,14,17,20,23,27,28,29,30,33,43,46,49,52,55,58,62,63,64,65,68,78,81,84,87,90,93,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,27,28,29,30,33,43,46,49,52,55,58,62,63,64,65,68,78,81,84,87,90,93,97,98,99,100,103),]
#for 12hr 9 live infections except P.sneebia + heatkilled
CountTable <- CountTable[,c(8,11,14,17,20,23,26,28,29,30,33,43,46,49,52,55,58,61,63,64,65,68,78,81,84,87,90,93,96,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,26,28,29,30,33,43,46,49,52,55,58,61,63,64,65,68,78,81,84,87,90,93,96,98,99,100,103),]

library("DESeq"); library(genefilter)
design$treatment <- factor(design$treatment)
attach(design)
cds.all = newCountDataSet(CountTable, design)
cds.all = estimateSizeFactors(cds.all)
cdsBlind = estimateDispersions(cds.all, method="blind", sharingMode="fit-only") #if replicates are not available, estimate across conditions.
vsd= varianceStabilizingTransformation(cdsBlind)


#PCA and then pull out the genes loading each PC
#The following method works similar to the one mentioned in https://www.biostars.org/p/13011/
rv = rowVars(exprs(vsd))
select = order(rv, decreasing=TRUE)[seq_len(500)]
#select = order(rv, decreasing=TRUE) #getting everything, not just the top 500
pca = prcomp(t(exprs(vsd)[select,]))
fac = factor(apply(pData(vsd)[, c("treatment"), drop = FALSE], 1, paste, collapse = " : "))

#proportion of variance
variance = pca$sdev^2 / sum(pca$sdev^2); variance = round(variance, 3) * 100
cat("Variance 1: ", variance[1], "Variance 2: ", variance[2])

#Get PC1 and PC2 genes and get the gene name
pca.rotation <- as.matrix(pca$rotation[,1:2]) #list of genes with PC1 and PC2
pca.rotation_ordered = pca.rotation[order(pca.rotation[,1], decreasing = T),]
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T) #Re-read full list of genes (subset of gtf file)
pca.rotation_ordered_genes = rownames(pca.rotation_ordered)
output.mx <- matrix(NA, ncol=1, nrow=length(pca.rotation_ordered_genes)); colnames(output.mx) <-c("gene_name")
for (i in 1:length(pca.rotation_ordered_genes)){
    gene.to.compare = as.character(pca.rotation_ordered_genes[i])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
pca.rotation_ordered = cbind(output.mx, pca.rotation_ordered)

#write.table(pca.rotation, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt", quote=F, row.names = T, col.names = T)
write.table(pca.rotation_ordered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/Answering_Nicolas_questions/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt", quote=F, row.names = T, col.names = T)
write.table(pca.rotation_ordered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/Answering_Nicolas_questions/genes.loaded.on.PC1.and.PC2.for.live.only.12hr.txt", quote=F, row.names = T, col.names = T)
write.table(pca.rotation_ordered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/Answering_Nicolas_questions/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.without.Staph.12hr.txt", quote=F, row.names = T, col.names = T)
write.table(pca.rotation_ordered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/Answering_Nicolas_questions/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.without.sneebia.12hr.txt", quote=F, row.names = T, col.names = T)


#Draw a PCA plot
pca.data = as.data.frame(pca$x)
X.ggplot =  pca.data[,1]; Y.ggplot =  pca.data[,2]; data = as.matrix(cbind(X.ggplot,Y.ggplot))
rownames(data) = rownames(pca.data)
treat = as.character(design$treatment)
treat=sub("S.mar.DB11", "S.marcescens Db11",treat); treat=sub("S.mar.type", "S.marcescens Type",treat)
treat=sub("E.faecalis.heat", "E.faecalis heatkilled",treat); treat=sub("P.rettgeri.heat", "P.rettgeri heatkilled",treat)
treat=sub("P. entomophila", "P.entomophila",treat)
design = cbind(treat,design[,2:3]); colnames(design) = c("treatment","time","rep")
PCs.with.design = cbind(data, design)
labels=rownames(PCs.with.design)

library(ggplot2)
cols = c("M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "S.aureus" = "orangered3","P.sneebia" = "blue",
         "S.marcescens Db11" = "navy", "P.entomophila" = "purple3", "E.faecalis heatkilled" = "dark salmon","P.rettgeri heatkilled" = "steelblue1")

#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_12hr_live_only_PCA_Apr_2016.pdf", height=10, width=13)
ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=8, aes(color = factor(treatment))) +
    scale_color_manual(values = cols, breaks=c("M.luteus", "E.coli", "S.marcescens Type","Ecc15","P.rettgeri","P.rettgeri heatkilled","E.faecalis","E.faecalis heatkilled","S.aureus","P.sneebia","S.marcescens Db11","P.entomophila"), name="Infection Conditions") + xlab(paste("PC1 (", variance[1], "%)",sep="")) + ylab(paste("PC2 (", variance[2], "%)",sep="")) + theme_bw() + theme(axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=20), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
#dev.off()


