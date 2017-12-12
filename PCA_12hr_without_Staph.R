#PCA analysis without S.aureus (Staph)
#October 16th, 2016
#Joo Hyun Im (ji72)

##Drawing a PCA with 12hr samples on laptop
rm(list=ls(all=TRUE)) #delete any previous entry

setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)

library("DESeq2")
#Only choosing 12hr samples (both live and heatkilled, but not clean prick, unchallenged, and S.aureus)
CountTable <- CountTable[,c(8,11,14,17,20,23,27,28,29,30,33,43,46,49,52,55,58,62,63,64,65,68,78,81,84,87,90,93,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,27,28,29,30,33,43,46,49,52,55,58,62,63,64,65,68,78,81,84,87,90,93,97,98,99,100,103),]
attach(design)

cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design, design = ~ treatment)
cds = estimateSizeFactors(cds)
cdsB = estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cdsB)

#pdf(file="Rplot.pdf")
#plotPCA(vsd, intgroup=c("treatment", "rep"))
#dev.off()

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
    
    cat("PC1 %: ", variance[1], ", PC2 %: ",variance[2])
    pca.data = as.data.frame(pca$x)
    X.ggplot =  pca.data[,1]; Y.ggplot =  pca.data[,2]; data = as.matrix(cbind(X.ggplot,Y.ggplot))
    rownames(data) = rownames(pca.data)
    cat("Variance 1: ", variance[1], "Variance 2: ", variance[2])
    write.table(data, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_live_only_PCA_without_Staph_Oct2016_only_PC1_and_PC2.txt", quote=F, col.names=T,row.names = T)
}
plotPCAWithSampleNames(vsd, intgroup=c("treatment"), ntop=500)
PCs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_live_only_PCA_without_Staph_Oct2016_only_PC1_and_PC2.txt", header= T) #A table of PC1 and PC2

#Change the treatment names to something more audience-friendly
treat = as.character(design$treatment)
treat=sub("S.mar.DB11", "S.marcescens Db11",treat); treat=sub("S.mar.type", "S.marcescens Type",treat)
treat=sub("E.faecalis.heat", "E.faecalis heatkilled",treat); treat=sub("P.rettgeri.heat", "P.rettgeri heatkilled",treat)
treat=sub("P. entomophila", "P.entomophila",treat)
design = cbind(treat,design[,2:3]); colnames(design) = c("treatment","time","rep")
#attach(design)
PCs.with.design = cbind(PCs, design)
labels=rownames(PCs.with.design)


library("ggplot2")
cols = c("M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "P.sneebia" = "blue", "S.marcescens Db11" = "navy", "P.entomophila" = "purple3", "E.faecalis heatkilled" = "dark salmon","P.rettgeri heatkilled" = "steelblue1")

pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq2_12hr_live_only_PCA_without_Staph_Oct2016.pdf", height=13, width=13)
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=8, aes(color = factor(treatment))) +
    scale_color_manual(values = cols, breaks=c("M.luteus", "E.coli", "S.marcescens Type","Ecc15","P.rettgeri","P.rettgeri heatkilled","E.faecalis","E.faecalis heatkilled","P.sneebia","S.marcescens Db11","P.entomophila"), name="Infection Conditions") + xlab("PC1 (42.0%)") + ylab("PC2 (16.7%)") + theme_bw() + theme(axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=20), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())

ggplot.PC
dev.off()
