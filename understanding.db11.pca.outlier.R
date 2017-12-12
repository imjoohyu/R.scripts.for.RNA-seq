#Understanding DB11 rep1: is it an outlier because of a few genes or is the entire replicate strange?
#April 10th, 2016
#Joo Hyun Im (ji72)

data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples.txt")
data.new = data[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)] #only 12hr samples

#group = c("M.luteus","E.coli","S.mar.type","Ecc15","P.rett","E.fae","Staph","P.sneebia","S.mar.DB11","PE","M.luteus","E.coli","S.mar.type","Ecc15","P.rett","E.fae","Staph","P.sneebia","PE","M.luteus","E.coli","S.mar.type","Ecc15","P.rett","E.fae","Staph","P.sneebia","S.mar.DB11","PE")
#########################################################
#Draw a correlation plot

# install.packages("corrplot")
library(corrplot)

#full data
cor.mat.data = cor(data, data)
corrplot(cor.mat.data, type="upper", tl.col = "black", tl.srt = 45, tl.cex=0.4)

#just Db11
cor.mat.db11 = cor(data[,c(28,63,98)],data[,c(28,63,98)])
corrplot(cor.mat.db11, type="upper", tl.col = "black", tl.srt = 45, tl.cex=0.4)

#data in raw number
cor.mat.raw = cor (data.new, data.new, method ="pearson")
corrplot(cor.mat.raw, type="upper", tl.col = "black", tl.srt = 45)

#data in log10
data.new.log10 = log10(data.new)
is.na(data.new.log10) = do.call(cbind,lapply(data.new.log10, is.infinite)) #turn all of the inf into NA
data.new.log10.na.omit = na.omit(data.new.log10)
cor.mat = cor(data.new.log10.na.omit, data.new.log10.na.omit, method="pearson")
corrplot(cor.mat, type="upper", order="hclust", tl.col="black", tl.srt=45)

#Checking correlation using r^2
FC.lm = lm(data.new[,9] ~ data.new[,21])
summary(FC.lm)$r.squared
FC.lm = lm(data.new[,9] ~ data.new[,33])
summary(FC.lm)$r.squared
FC.lm = lm(data.new[,21] ~ data.new[,33])
summary(FC.lm)$r.squared

#What if I used the filtered dataset to draw PCA?
#I created an MDS plot after filtering out for low-count genes and the MDS plot looks the same. ID_28 still sticks out (see Evernote note called "What's gong on with Db11?")
#Still the same

#########################################################
#Is this pattern driven by a few genes?
#1. Find the genes that are outlier.
data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples.txt")
data.db11 = data[,c(28,63,98)]

#Q. Does the log10 transformation improve the R^2 values?
data.db11.log10 = log10(data.db11) #data.db11.log10 = log10(data.db11.norm) is the same as unnormalized
is.na(data.db11.log10) = do.call(cbind,lapply(data.db11.log10, is.infinite)) #turn all of the inf into NA
data.db11.log10.na.omit = na.omit(data.db11.log10)
summary(lm(data.db11.log10.na.omit[,1] ~ data.db11.log10.na.omit[,2]))$r.squared #28 & 63 = 0.942
summary(lm(data.db11.log10.na.omit[,1] ~ data.db11.log10.na.omit[,3]))$r.squared #28 & 98 = 0.933
summary(lm(data.db11.log10.na.omit[,2] ~ data.db11.log10.na.omit[,3]))$r.squared #63 & 98 = 0.927
#A. Yes, log10 transformation improves the R^2 values.

#Q. Does normalizing the replicates by library size change the pattern?
data.db11.norm = data.frame(cbind((data.db11[,1]/1.10), data.db11[,2], data.db11[,3])) #normalize it by size
colnames(data.db11.norm) = colnames(data.db11)
rownames(data.db11.norm) = rownames(data.db11)
plot(log10(data.db11), main="Unnormalized")
#plot(log10(data.db11.norm), main= "Normalized") 
#A. No, not much different

#Q. Is one point skewing the R^2 value?
data.db11.ID28.sorted = data.db11[order(-data.db11$ID_28),] #on raw count
data.db11.ID28.sorted.remove.max = data.db11.ID28.sorted[2:17558,] #remove the biggest count that skews the correlation
summary(lm(data.db11.ID28.sorted.remove.max[,1] ~data.db11.ID28.sorted.remove.max[,2]))$r.squared #28 & 63 = 0.944 #his improved tremendously from 0.63
summary(lm(data.db11.ID28.sorted.remove.max[,1] ~data.db11.ID28.sorted.remove.max[,3]))$r.squared #28 & 98 = 0.915
summary(lm(data.db11.ID28.sorted.remove.max[,2] ~data.db11.ID28.sorted.remove.max[,3]))$r.squared #63 & 98 = 0.860 #This improved tremendously from 0.58
cor.mat.max.removed = cor(data.db11.ID28.sorted.remove.max,data.db11.ID28.sorted.remove.max)
corrplot(cor.mat.max.removed, type="upper", tl.col = "black", tl.srt = 45)
#A. Getting rid of the gene with the highest number of count improved the R^2 values.

#Q. Does ID_63 correlate well with other samples after having the max gene (FBgn0013686) taken out?
data.new.max.removed = data.new[order(-data.new$ID_8),]
data.new.max.removed = data.new.max.removed[2:17558,] #removing FBgn0013686
cor.mat.all.max.removed = cor(data.new.max.removed, data.new.max.removed)
corrplot(cor.mat.all.max.removed, type="upper", tl.col = "black", tl.srt = 45, tl.cex=0.4, cl.cex = 0.5)
#A. It looks more normal.

plot(log10(data.db11.ID28.sorted.remove.max), main="Unnormalized, max point removed")

#Q. Does taking out of the genes that have the highest variance among replicates imrpove the correlation?
#Now pick genes that are deviating in ID_28 compared to other replicates
data.db11.var = apply(data.db11, 1, var)
data.db11.with.var = cbind(data.db11,data.db11.var)
data.db11.with.var.ordered = data.db11.with.var[order(-data.db11.with.var$data.db11.var),]
plot(log10(data.db11.with.var.ordered[6:17558,1:3]), main="Db11, without the top most variable genes") #Get rid of the top 5
cor.mat.all.max.5.removed = cor(data.db11.with.var.ordered[6:17558,1:3],data.db11.with.var.ordered[6:17558,1:3])
corrplot(cor.mat.all.max.5.removed, type="upper", tl.col = "black", tl.srt = 45, tl.cex=0.4, cl.cex = 0.5)

#2. Remove them from the whole data (for now) and re-run the PCA. Does PCA change?
count.table = CountTable[,c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99)] 
#remove = c("FBgn0013686","FBgn0000276","FBgn0000277")
remove = row.names(data.db11.with.var.ordered[1:20,])
count.table2 = count.table[!rownames(count.table) %in% remove, ] #Remove FBgn0013686, FBgn0000276 CecA1, and FBgn0000277 CecA2
group = c("Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE") 
colors = c("lightskyblue","darkseagreen2","coral1","lightpink","dodgerblue2","blue","purple3","lightseagreen","navy","orangered3")
group = factor(group); points= c(0,1,2,3,4,5,6,7,8,9)

d1 = DGEList(counts = count.table, group=group); d2 = DGEList(counts = count.table2, group=group); 
d1 = calcNormFactors(d1); d2 = calcNormFactors(d2)
#pdf(file="/home/ji72/MDS_12hr_only_all_missing_CecA1_CecA2_Apr13_2016.pdf", width=15, height=10)
#pdf(file="/Users/JooHyun/Desktop/MDS_12hr_only_all_missing_CecA1_CecA2_Apr13_2016.pdf", width=15, height=10)
pdf(file="/Users/JooHyun/Desktop/MDS_12hr_only_all_missing_top20_variable_genes_Apr13_2016.pdf", width=15, height=10)
par(mfrow=c(1,2))
plotMDS(d1, col=colors[group], pch=points[group], main="12hr live", cex=1.5)
legend("topright", legend=levels(group), pch=points, col=colors, ncol=2, pt.cex=1)
plotMDS(d2, col=colors[group], pch=points[group], main="12hr live minus 20 variable genes", cex=1.5)
legend("topright", legend=levels(group), pch=points, col=colors, ncol=2, pt.cex=1)
dev.off()
#Take out more genes (like top 20 variable genes). Taking out three/twenty genes didn't really change the PCA shape.

plot(log10(count.table2[,c(9,19,29)]), main = "Db11, without the top 20 most variable genes")

#Q. Does the ID_28 (Db11 rep1) correlate better with P. sneebia or S.marcescens Type than other Db11 replicates?
library(corrplot)
data.db11.sneb.smar.log10 = log10(data[,c(28,63,98,27,62,97,14,49,84,20,55,90 )]) #Db11, Psneeb, Smartype, Prett
is.na(data.db11.sneb.smar.log10) = do.call(cbind,lapply(data.db11.sneb.smar.log10, is.infinite)) #turn all of the inf into NA
data.db11.sneb.smar.log10.na.omit = na.omit(data.db11.sneb.smar.log10)
cor.mat.data.db11.sneb.smar = cor(data.db11.sneb.smar.log10.na.omit,data.db11.sneb.smar.log10.na.omit)
corrplot(cor.mat.data.db11.sneb.smar, method="number",type="upper", tl.col = "black", tl.srt = 45)
#Answer: ID_28 is highly correlated with both its replicates (0.97) and other species (0.97-0.98) at least at the log10 level.
#I don't think the outlier pattern of ID_28 is due to a few genes.

#Q.Find out which genes contribute most to the PCA PC1
#Correlate PC1 and gene expression (PC1 ~ gene1 exp across IDs) = number of genes, 
#order them and find the highest as the genes that contribute most to PCA (lowest p-val = most association with PC1)
data.new = data[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)] #only 12hr samples
PC = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_only_PC1_and_PC2_April_2016.txt")
#cor.between.PC1.and.data = apply(c(data.new,PC[,1]),1,cor)
data.new.log10 = log10(data.new)
is.na(data.new.log10) = do.call(cbind,lapply(data.new.log10, is.infinite)) #turn all of the inf into NA
data.new.log10.na.omit = na.omit(data.new.log10)

pval.table = matrix(data = NA, nrow=dim(data.new.log10.na.omit)[1], ncol=2)
for (x in 1:dim(data.new.log10.na.omit)[1]){
    lm.pca.pvalue = summary(lm(t(data.new.log10.na.omit[x,1:36]) ~ PC[1:36,1]))$coef[,4]
    pval.table[x,1] = rownames(data.new.log10.na.omit[x,])
    pval.table[x,2] = lm.pca.pvalue[2]
}

pval.table.sorted = pval.table[order(pval.table[,2]),]
pval.table.sorted.sig = pval.table[which(pval.table[,2] < 0.05),]
write.table(pval.table.sorted.sig[,1], file="/Users/JooHyun/Desktop/genes.enriched.in.12hr.PCA.txt", quote=F, row.names = F)
#A. Genes (with pval <0.05) are enriched in anatomical structure development, positive regulation of transcription from RNA polymerase II promoter, etc


#########################################################
#PCA: Getting rid of an outlier replicate
#What if I exclude ID_63? How about excluding ID_28?
rm(list=ls(all=TRUE)) #delete any previous entry
library("edgeR")
CountTable =read.table("mega_RNA-seq_count_of_all_samples.txt")
count.table1 = CountTable[,c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99)] #all 30 samples
count.table2 = CountTable[,c(8,11,14,17,20,23,26,27,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99)] #ID_28 missing
count.table3 = CountTable[,c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,64,78,81,84,87,90,93,96,97,98,99)] #ID_63 missing
group1 = c("Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE") 
group2 = c("Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","PE",
            "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
            "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE") #ID_28 missing
group3 = c("Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","PE",
           "Mluteus","Ecoli","Smar.type","Efaecalis","Ecc15","P.rett","Staph","P.sneebia","SmarDB11","PE") #ID_63 missing
group1 = factor(group1); group2 = factor(group2); group3 = factor(group3)

colors = c("lightskyblue","darkseagreen2","coral1","lightpink","dodgerblue2","blue","purple3","lightseagreen","navy","orangered3")
points= c(0,1,2,3,4,5,6,7,8,9)
d1 = DGEList(counts = count.table1, group=group1); d2 = DGEList(counts = count.table2, group=group2); d3 = DGEList(counts = count.table3, group=group3)
d1 = calcNormFactors(d1); d2 = calcNormFactors(d2); d3 = calcNormFactors(d3)

#pdf(file="/home/ji72/MDS_without_ID_63.pdf")
pdf(file="/home/ji72/MDS_12hr_only_all_rep1_missing_rep2_missing.pdf", width=17, height=10)
par(mfrow=c(1,3))
plotMDS(d1, col=colors[group1], pch=points[group1], main="12hr live", cex=1.5)
legend("topright", legend=levels(group1), pch=points, col=colors, ncol=2, pt.cex=1)
plotMDS(d2, col=colors[group2], pch=points[group2], main="12hr live missing Db11 rep1 (ID_28)", cex=1.5)
legend("topright", legend=levels(group2), pch=points, col=colors, ncol=2, pt.cex=1)
plotMDS(d3, col=colors[group3], pch=points[group3], main="12hr live missing Db11 rep2 (ID_63)", cex=1.5)
legend("topright", legend=levels(group3), pch=points, col=colors, ncol=2, pt.cex=1)
dev.off()

pdf(file="/home/ji72/MDS_12hr_only_all_rep1_missing_rep2_missing_label.pdf", width=17, height=10)
par(mfrow=c(1,3))
plotMDS(d1); plotMDS(d2); plotMDS(d3)
dev.off()

#Answer: Removing ID_28 or ID_63 doesn't really change the rest of the PCA


#Q. What are the genes that drive Db11 replicates apart?

data.db11 = data[,c(28,63,98)]

rv = rowVars(assay(x))
select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(x)[select,]))

# proportion of variance
variance = pca$sdev^2 / sum(pca$sdev^2)
variance = round(variance, 3) * 100

# sample names
names = colnames(x)

# factor of groups
fac.time = factor(apply(as.data.frame(colData(x)[, intgroup[1], drop=FALSE]), 1, paste, collapse=" : "))
print(fac.time)
fac.treatment = factor(apply(as.data.frame(colData(x)[, intgroup[2], drop=FALSE]), 1, paste, collapse=" : "))
