#edgeR: comparing S.mar.type vs Db11 and E.fae.live vs E.fae.hk
#November 5, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
CountTable.s = CountTable[,c(3:104)]
library("edgeR")

S.mar.type.12hr = CountTable.s[,c(13,47,81)]; S.mar.Db11.12hr = CountTable.s[,c(27,61,95)]
E.fae.live.12hr = CountTable.s[,c(16,50,84)]; E.fae.live.36hr = CountTable.s[,c(17,51,85)]; E.fae.live.5.5d = CountTable.s[,c(18,52,86)]
E.fae.heatkilled.12hr = CountTable.s[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable.s[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable.s[,c(31,65,99)]

full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T) #Re-read full list of genes (subset of gtf file)

#1. S.marcescens Type vs S.marcescens Db11
count.table = cbind(S.mar.type.12hr,S.mar.Db11.12hr); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("Type",3), rep("Db11",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("Type",3), rep("Db11",3))); design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
et = exactTest(de, pair=c("Type", "Db11")); tt = topTags(et,n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #9 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in Db11
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #9 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #1 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/S.mar.type.12h.vs.S.mar.Db11.12h.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/S.mar.type.12h.vs.S.mar.Db11.12h.down.txt", quote=F, row.names = F, col.names = T)
#Results: 9 DEGs that are upregulated in Db11 infection; 1 DEG that are downregulated in Db11 infection

#2. E.faecalis live 12hr vs E.faecalis heatkilled 12hr
count.table = cbind(E.fae.live.12hr,E.fae.heatkilled.12hr); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("E.fae.live",3), rep("E.fae.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("E.fae.live",3), rep("E.fae.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
et = exactTest(de, pair=c("E.fae.hk", "E.fae.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #182 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in E.faecalis live
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
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #95 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.12h.vs.E.fae.hk.12h.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.12h.vs.E.fae.hk.12h.down.txt", quote=F, row.names = F, col.names = T)
#Results: 88 DEGs that are upregulated in live bacterial infection; 95 DEGs that are downregulated in live bacterial infection

#3. E.faecalis live 36hr vs E.faecalis heatkilled 36hr
count.table = cbind(E.fae.live.36hr,E.fae.heatkilled.36hr); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("E.fae.live",3), rep("E.fae.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("E.fae.live",3), rep("E.fae.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
et = exactTest(de, pair=c("E.fae.hk", "E.fae.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #84 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in E.faecalis live
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #63 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #21 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.36h.vs.E.fae.hk.36h.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.36h.vs.E.fae.hk.36h.down.txt", quote=F, row.names = F, col.names = T)
#Results: 63 DEGs that are upregulated in live bacterial infection; 21 DEGs that are downregulated in live bacterial infection

#4. E.faecalis live 5.5d vs E.faecalis heatkilled 5.5d
count.table = cbind(E.fae.live.5.5d,E.fae.heatkilled.5.5d); rownames(count.table) = CountTable[,1]
design = data.frame(row.names = colnames(count.table), condition = c(rep("E.fae.live",3), rep("E.fae.hk",3)))
d = DGEList(counts = count.table, group=design$condition); d = calcNormFactors(d) #Normalization
group = factor(c(rep("E.fae.live",3), rep("E.fae.hk",3)))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
#Find any genes that differ between any of the treatment conditions 12hr and 36hr
et = exactTest(de, pair=c("E.fae.hk", "E.fae.live")); tt = topTags(et, n=nrow(de)) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05]; plotSmear(de, de.tags=deg, cex=1) # MA plot
sig.genes = tt$table[which(tt$table$FDR < 0.05),] #46 genes
sig.genes.ordered.by.logFC = sig.genes[order(as.numeric(sig.genes$logFC), decreasing=T),] #Positive log2FC = upregulation in E.faecalis live
sig.genes.ordered.by.logFC = cbind(rownames(sig.genes.ordered.by.logFC),sig.genes.ordered.by.logFC) 
list.of.unknown.genes = sig.genes.ordered.by.logFC; colnames(list.of.unknown.genes) <- c("gene_id")
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1]); colnames(output.mx) <-c("gene_name")
for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}
with.name = cbind(sig.genes.ordered.by.logFC[,1], output.mx, sig.genes.ordered.by.logFC[,c(2:5)]); colnames(with.name) = c("gene_id","gene_name","logFC","logCPM","P-value","FDR")
sig.genes.ordered.by.logFC.upr = with.name[which(with.name$logFC > 0),] #46 genes
sig.genes.ordered.by.logFC.downr = with.name[which(with.name$logFC < 0),] #0 genes
write.table(sig.genes.ordered.by.logFC.upr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.5.5d.vs.E.fae.hk.5.5d.up.txt", quote=F, row.names = F, col.names = T)
write.table(sig.genes.ordered.by.logFC.downr, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/E.fae.live.5.5d.vs.E.fae.hk.5.5d.down.txt", quote=F, row.names = F, col.names = T)
#Results: 46 DEGs that are upregulated in live bacterial infection; 0 DEGs that are downregulated in live bacterial infection
