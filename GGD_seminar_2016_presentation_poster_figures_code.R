#Making a PCA plot featuring 12hr time points of live bacteria on CBSU
#February 2nd, 2016
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)

library("DESeq2")
#Only choosing 12hr samples (All live bacteria, 12hr unmolested, 12hr clean prick // Heat-killed samples were excluded)
CountTable <- CountTable[,c(2,5,8,11,14,17,20,23,26,27,28,29,37,40,43,46,49,52,55,58,61,62,63,64,72,75,78,81,84,87,90,93,96,97,98,99)]
design <- design[c(2,5,8,11,14,17,20,23,26,27,28,29,37,40,43,46,49,52,55,58,61,62,63,64,72,75,78,81,84,87,90,93,96,97,98,99),]
attach(design)

#Only choosing 12hr samples (All live bacteria // 12hr unmolested, 12hr clean prick and heat-killed samples were excluded)
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99),]
attach(design)

cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design, design = ~ treatment)
cds = estimateSizeFactors(cds)
cdsB = estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cdsB)
library(RColorBrewer)
library(genefilter)
library(lattice)
plotPCAWithSampleNames = function(x, intgroup=c("treatment"), ntop=500) {
    # pca
    rv = rowVars(assay(x))
    select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select,]))
    
    # proportion of variance
    variance = pca$sdev^2 / sum(pca$sdev^2)
    variance = round(variance, 3) * 100
    
    # sample names
    names = colnames(x)
    
    # factor of groups
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
    
    # colors
    if( nlevels(fac) >= 10 )
        colors = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"), brewer.pal(12,"Set3"), brewer.pal(7,"Dark2"))
    
    # plot
    xyplot(
        PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
        aspect = "fill",
        col = colors,
        xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
        ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
        panel = function(x, y, ...) {
            panel.xyplot(x, y, ...);
            ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
        },
        main = draw.key(
            key = list(
                rect = list(col = colors),
                text = list(levels(fac)),
                rep = FALSE
            )
        )
    )
    
    cat("PC1 %: ", variance[1], ", PC2 %: ",variance[2], ", PC3 %: ",variance[3])
    pca.data = as.data.frame(pca$x)
    write.table(pca.data[,c(1:3)], file="/home/ji72/mega_RNA-seq_DESeq2_12hr_live_only_PC1-3.txt", quote=F, col.names=T,row.names = T)
}
plotPCAWithSampleNames(vsd, intgroup=c("treatment"), ntop=500)


#On local machine
#PCs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_live_only_PC1_and_PC2.txt", header= T) #A table of PC1 and PC2
PCs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_live_only_PC1-3.txt", header= T) #A table of PC1, PC2, and PC3
design <- read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)
design <- design[c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99),]
#design <- design[c(2,5,8,11,14,17,20,23,26,27,28,29,37,40,43,46,49,52,55,58,61,62,63,64,72,75,78,81,84,87,90,93,96,97,98,99),]



#Change the treatment names to something more audience-friendly
treat = as.character(design$treatment)
treat = sub("Unmolested", "Uninjured",treat); treat=sub("Clean.prick", "Sterile wounded",treat)
treat=sub("S.mar.DB11", "S.marcescens Db11",treat); treat=sub("S.mar.type", "S.marcescens Type",treat)
design = cbind(treat,design[,2:3]); colnames(design) = c("treatment","time","rep")
attach(design)
PCs.with.design = cbind(PCs, design)
labels=rownames(PCs.with.design)

library("ggplot2")

#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_12hr_live_only_PCA_Feb_2016.pdf", height=10, width=13) 

#2/1: change colors of the dots 
cols = c("Uninjured" = "azure3", "Sterile wounded" = "azure4", "M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "S.aureus" = "orangered3","P.sneebia" = "blue","S.marcescens Db11" = "navy", "P.entomophila" = "purple3")
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=7, aes(color = factor(treatment))) + scale_color_manual(values=cols, breaks=c("Uninjured", "Sterile wounded",
       "M.luteus", "E.coli", "S.marcescens Type","Ecc15","P.rettgeri","E.faecalis","S.aureus","P.sneebia","S.marcescens Db11","P.entomophila"), name="Infection Condition") + xlab("PC1 (40.0%)") + ylab("PC2 (21.9%)") +
    theme_bw() + theme(legend.justification=c(0,0), legend.position=c(0,0), axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=20), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) #+ xlim(-20,20) + ylim(-20,20)
ggplot.PC


#2/16: making a 3D PCA graph
cols = c("M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "coral1", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "lightskyblue", "S.aureus" = "orangered3","P.sneebia" = "blue","S.marcescens Db11" = "navy", "P.entomophila" = "purple3")  # had to manually swap between Ecc15 and E.faecalis colors because they were swapped automatically.
library(rgl)
attach(PCs.with.design)
plot3d(PC1,PC2,PC3,xlab = "PC1 (39.3%)", ylab="PC2 (23.2%)", zlab="PC3 (8.0%)", col=cols, size = 10)
#with(PCs.with.design, plot3d(PC1, PC2, PC3)
with(PCs.with.design, text3d(PC1, PC2, PC3, treatment))

#==============================================================================================================
#2. Count the number of DEGs per infection-time condition (Unchallenged vs Infected)
#read the file that has a number of DEGs per infection-time condition
num.degs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/num.of.DEGs.Nov.2015.txt", header=F)
#Number from list_of_DEGs folder created by GO.analysis.pulling.out.sig.genes.
num.degs.down = num.degs[1:31,]; num.degs.up = num.degs[32:62,]
par(mfrow=c(2,3)) 
cols2 = c("Sterile wounded" = "azure4", "M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "S.aureus" = "orangered3","P.sneebia" = "blue","S.marcescens Db11" = "navy", "P.entomophila" = "purple3")
barplot(num.degs.up[c(1:6,8,10:13),2], col=cols2, ylim=c(0,700), xlab="12hr",cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.up[c(14:19,21),2], col=cols2, ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.up[c(23:28,30),2], col=cols2, ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot
barplot(num.degs.down[c(1:6,8,10:13),2], col=cols2, ylim=c(0,700), xlab="12hr", cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.down[c(14:19,21),2], col=cols2, ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.down[c(23:28,30),2], col=cols2, ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot

#Just get the legend
par(mfrow=c(1,1)) 
barplot(c(1:11), col=c("azure4","lightpink","darkseagreen2","lightseagreen","lightskyblue","dodgerblue2","coral1","orangered3","blue","navy","purple3"), 
        legend = c("Sterile wounded","M.luteus","E.coli","S.marcescens Type","Ecc15","P.rettgeri",
                   "E.faecalis", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila"))

#Clean-prick vs Infected
#read the file that has a number of DEGs per infection-time condition
num.degs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_comp_with_cleanprick_Dec_2015/list_of_DEGs/num.of.DEGs.clean.prick.Dec.2015.txt", header=F)
#Number from list_of_DEGs folder created by GO.analysis.pulling.out.sig.genes.
num.degs.down = num.degs[1:28,]; num.degs.up = num.degs[29:56,]
par(mfrow=c(2,3)) 
cols2 = c("M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
          "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "S.aureus" = "orangered3","P.sneebia" = "blue","S.marcescens Db11" = "navy", "P.entomophila" = "purple3")
barplot(num.degs.up[c(1:5,7,9:12),2], col=cols2, ylim=c(0,700), xlab="12hr",cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.up[c(13:17,19),2], col=cols2, ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.up[c(21:25,27),2], col=cols2, ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot
barplot(num.degs.down[c(1:5,7,9:12),2], col=cols2, ylim=c(0,700), xlab="12hr", cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.down[c(13:17,19),2], col=cols2, ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.down[c(21:25,27),2], col=cols2, ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot

#Just get the legend
par(mfrow=c(1,1)) 
barplot(c(1:11), col=c("azure4","lightpink","darkseagreen2","lightseagreen","lightskyblue","dodgerblue2","coral1","orangered3","blue","navy","purple3"), 
        legend = c("Sterile wounded","M.luteus","E.coli","S.marcescens Type","Ecc15","P.rettgeri",
                   "E.faecalis", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila"))

#3. Staph-unique response compared to M.luteus infeciton at 12hr
require(ggplot);require(RColorBrewer)
bacteria = c("S.aureus","S.aureus"); bacteria2 = c("M.luteus","M.luteus")
names = c("Unique","All DEGs"); names2 = c("Unique","All DEGs")
numbers = c(29,71); numbers2 = c(5.9,94.1) #put the overlap number first
df <- data.frame(bacteria,names,numbers); df2 <- data.frame(bacteria2,names2,numbers2)


ggplot(df, aes(x = bacteria, width=10, height=7)) +geom_bar(aes(weight=numbers, fill = rev(names), position = 'stack')) + 
    scale_y_continuous("Percentage of DEGs that are unique to pathogen") + scale_x_discrete(NULL) + scale_fill_manual(values = c("orangered3","ivory")) +
    theme(axis.title = element_text(face = "bold", size = 30), legend.title = element_text(size = 30), legend.text = element_text(size=30), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme_bw(base_size = 35)
ggplot(df2, aes(x = bacteria2, width=10, height=7)) +geom_bar(aes(weight=numbers2, fill = rev(names2), position = 'stack')) + 
    scale_y_continuous("Percentage of DEGs that are unique to pathogen") + scale_x_discrete(NULL) + scale_fill_manual(values = c("lightpink","ivory")) +
    theme(axis.title = element_text(face = "bold", size = 30), legend.title = element_text(size = 30), legend.text = element_text(size=30), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme_bw(base_size = 35)




