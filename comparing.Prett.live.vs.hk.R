#edgeR: comparing P.rett.live vs P.rett.hk
#December 7, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
CountTable.s = CountTable[,c(3:104)]
library("edgeR")

P.rett.live.12hr = CountTable.s[,c(19,53,87)]; P.rett.live.36hr = CountTable.s[,c(20,54,88)]; P.rett.live.5.5d = CountTable.s[,c(21,55,89)]
P.rett.heatkilled.12hr = CountTable.s[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable.s[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable.s[,c(34,68,102)]

full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T) #Re-read full list of genes (subset of gtf file)

#1. P.rettgeri live 12hr vs P.rettgeri heatkilled 12hr
count.table = cbind(P.rett.live.12hr,P.rett.heatkilled.12hr); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("P.rett.live",3), rep("P.rett.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("P.rett.live",3), rep("P.rett.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
et = exactTest(de, pair=c("P.rett.hk", "P.rett.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #344 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in P.rettcalis live
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #117 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #227 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.12h.vs.P.rett.hk.12h.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.12h.vs.P.rett.hk.12h.down.txt", quote=F, row.names = F, col.names = T)
#Results: 117 DEGs that are upregulated in live bacterial infection; 227 DEGs that are downregulated in live bacterial infection

#3. P.rettgeri live 36hr vs P.rettgeri heatkilled 36hr
count.table = cbind(P.rett.live.36hr,P.rett.heatkilled.36hr); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("P.rett.live",3), rep("P.rett.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("P.rett.live",3), rep("P.rett.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
et = exactTest(de, pair=c("P.rett.hk", "P.rett.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #298 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in P.rettcalis live
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #152 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #146 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.36h.vs.P.rett.hk.36h.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.36h.vs.P.rett.hk.36h.down.txt", quote=F, row.names = F, col.names = T)
#Results: 152 DEGs that are upregulated in live bacterial infection; 146 DEGs that are downregulated in live bacterial infection

#4. P.rettgeri live 5.5d vs P.rettgeri heatkilled 5.5d
count.table = cbind(P.rett.live.5.5d,P.rett.heatkilled.5.5d); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("P.rett.live",3), rep("P.rett.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("P.rett.live",3), rep("P.rett.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
#Find any genes that differ between any of the treatment conditions 12hr and 36hr
et = exactTest(de, pair=c("P.rett.hk", "P.rett.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #108 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in P.rettcalis live
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #88 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #20 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.5.5d.vs.P.rett.hk.5.5d.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/P.rett.live.5.5d.vs.P.rett.hk.5.5d.down.txt", quote=F, row.names = F, col.names = T)
#Results: 88 DEGs that are upregulated in live bacterial infection; 20 DEGs that are downregulated in live bacterial infection
