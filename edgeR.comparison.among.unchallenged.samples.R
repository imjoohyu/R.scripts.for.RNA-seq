#edgeR: comparing unchallenged samples (12hr, 36hr, 5.5d)
#November 4, 2015
#Joo Hyun Im (ji72)

#*****Conclusion: We will exclude 0hr unchallenged sample and use 12hr, 36hr, and 5.5day unchallenged samples as 9 replicates for control (2/5/2016)

rm(list=ls(all=TRUE))

CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_with_gene_id_and_names.txt", header=T) #11558 x 107
library("edgeR")

#1. Comparing three time points of the unchallenged samples to see if they are radically different
CountTable.uc = CountTable[,c(4,5,6,39,40,41,74,75,76)] # total 9 samples
rownames(CountTable.uc) = CountTable[,1]
design = data.frame(row.names = colnames(CountTable.uc), condition = c("twelve","thirty.six","five.half","twelve","thirty.six","five.half","twelve","thirty.six","five.half"))
d = DGEList(counts = CountTable.uc, group=design$condition); keep = rowSums(cpm(d) > 1.2) >=3 #at least 1/3 of the samples
summary(keep); d = d[keep, , keep.lib.sizes=FALSE]; d = calcNormFactors(d); plotMDS(d) #Normalization
group = factor(c("twelve","thirty.six","five.half","twelve","thirty.six","five.half","twelve","thirty.six","five.half"))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
#Find any genes that differ between any of the treatment conditions 12hr, 36hr, or 5.5d
#"Technically, this procedure tests whether either of the contrats" 36hr-12hr or 5.5d-12hr are non-zero (ANOVA-like)
fit = glmFit(de, design.mx); lrt = glmLRT(fit, coef=2:3) #coef=2:3 because I need to exclude the intercept
tt = topTags(lrt) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR"); rownames(de.with.specific.info) = rownames(tt.ordered.by.gene.id)
#Results: there are only 10 DEGs among three time points of the unchallenged samples


#2. Comparing three time points of the unchallenged samples to see if they are radically different after removing two samples that didn't cluster together
#MDS plot showed that ID_38 and ID_72 are not close to each other. What if we remove them and do it again? (there are one 12hr and one 36hr unchallenged samples)
CountTable.uc.removed = CountTable[,c(4,5,6,39,41,75,76)]; rownames(CountTable.uc.removed) = CountTable[,1]
design = data.frame(row.names = colnames(CountTable.uc.removed), condition = c("twelve","thirty.six","five.half","twelve","five.half","thirty.six","five.half"))
d = DGEList(counts = CountTable.uc.removed, group=design$condition); keep = rowSums(cpm(d) > 1.2) >=2 #2 out of 7
summary(keep); d = d[keep, , keep.lib.sizes=FALSE]; d = calcNormFactors(d)#Normalization
group = factor(c("twelve","thirty.six","five.half","twelve","five.half","thirty.six","five.half"))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
#Find any genes that differ between any of the treatment conditions 12hr, 36hr, or 5.5d
#"Technically, this procedure tests whether either of the contrats" 36hr-12hr or 5.5d-12hr are non-zero (ANOVA-like)
fit = glmFit(de, design.mx); lrt = glmLRT(fit, coef=2:3) #coef=2:3 because I need to exclude the intercept
tt = topTags(lrt) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR"); rownames(de.with.specific.info) = rownames(tt.ordered.by.gene.id)
#Results: there are only 10 DEGs among three time points after removing the unchallenged samples that did not cluster with each other.
#Most (8 out of 10) genes overlapped. So there are still about 10 genes that are DE among unchallenged samples.
#In other words, there is no point in removing them.

#3. What if we look at UC-12hr, 36hr, and 5.5d with 0hr? Do all of the 12,36,5.5 cluster together and separate from 0hr?
CountTable.uc.all = CountTable[,c(3,4,5,6,38,39,40,41,73,74,75,76)] # total 12 samples
rownames(CountTable.uc.all) = CountTable[,1]
design = data.frame(row.names = colnames(CountTable.uc.all), condition = c("zero","twelve","thirty.six","five.half","zero", "twelve","thirty.six","five.half","zero","twelve","thirty.six","five.half"))
d = DGEList(counts = CountTable.uc.all, group=design$condition); keep = rowSums(cpm(d) > 1.2) >=4 #at least 1/3 of the samples
summary(keep); d = d[keep, , keep.lib.sizes=FALSE]; d = calcNormFactors(d); plotMDS(d) 
#Results: Here, two of the 0hr samples (ID_1, and ID_71) are not clustering with others as expected.
#However, one sample of 12hr unchallenged (ID_72) also does not cluster with the rest. This happened when just looking at UC-12hr,36hr,5.5d only.
#We could re-run the analyses without the ID_72. But without ID_72, there will be a few DEGs among unchallenged samples (regardless of time points).
#So that may not be useful.

#4. Are there any differences among 12hr and 36hr?
CountTable.uc.two = CountTable[,c(4,5,39,40,74,75)] # total 6 samples
design = data.frame(row.names = colnames(CountTable.uc.two), condition = c("twelve","thirty.six","twelve","thirty.six","twelve","thirty.six"))
d = DGEList(counts = CountTable.uc.two, group=design$condition); keep = rowSums(cpm(d) > 1.2) >=2 #2 out of 6
summary(keep); d = d[keep, , keep.lib.sizes=FALSE]; d = calcNormFactors(d)#Normalization
group = factor(c("twelve","thirty.six","twelve","thirty.six","twelve","thirty.six"))
design.mx = model.matrix(~group, data=d$samples); de = estimateDisp(d, design.mx) #Estimate dispersion
#Find any genes that differ between any of the treatment conditions 12hr and 36hr
et = exactTest(de); tt = topTags(et) #table, adjust.method, comparison, test
rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR"); rownames(de.with.specific.info) = rownames(tt.ordered.by.gene.id)
#Results: Yes, there will always be DEGs.

#5. Is there any pattern in the way these DEGs are expressed over time? (from #1)
CountTable.uc.degs.only = CountTable.uc[rownames(tt.ordered.by.gene.id),]
log2.CountTable.uc.degs.only = log2(CountTable.uc.degs.only)
twelve = apply(cbind(log2.CountTable.uc.degs.only[,1],log2.CountTable.uc.degs.only[,4],log2.CountTable.uc.degs.only[,7]),1,mean)
thirty.six = apply(cbind(log2.CountTable.uc.degs.only[,2],log2.CountTable.uc.degs.only[,5],log2.CountTable.uc.degs.only[,8]),1,mean)
five.half = apply(cbind(log2.CountTable.uc.degs.only[,3],log2.CountTable.uc.degs.only[,6],log2.CountTable.uc.degs.only[,9]),1,mean)
total = cbind(twelve, thirty.six, five.half)


clean.prick = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs_in_each_conditions/list_of_upregulated_DE_genes_in_clean.prick_name_only.txt", header=T)
write.table(unique(clean.prick), file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs_in_each_conditions/list_of_upregulated_DE_genes_in_clean.prick_name_only_uniq.txt", quote=F, row.names = F)
