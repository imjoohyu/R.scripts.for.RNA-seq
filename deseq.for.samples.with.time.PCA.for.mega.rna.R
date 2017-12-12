#DESeq Commands for Mega RNA-seq: PCA
#Date: September 2, 2015
#ji72

rm(list=ls(all=TRUE)) #delete any previous entry

#After reading all the samples into R
#mega.rnaseq.read.table.R

#All together
list.of.samples <- vector()
list.of.samples.without.comma <- vector()
for (i in 1:105){
    list.of.samples[i] <- paste("ID_",i,sep="")
    list.of.samples.without.comma[i] <- paste("ID_",i,",",sep="")
}

#Put all the dataset together
countTable <-cbind(ID_1, ID_2, ID_3, ID_4, ID_5, ID_6, ID_7, ID_8, ID_9, ID_10, ID_11, ID_12, ID_13, ID_14, ID_15, ID_16, 
                   ID_17, ID_18, ID_19, ID_20, ID_21, ID_22, ID_23, ID_24, ID_25, ID_26, ID_27, ID_28, ID_29, ID_30, ID_31,
                   ID_32, ID_33, ID_34, ID_35, ID_36, ID_37, ID_38, ID_39, ID_40, ID_41, ID_42, ID_43, ID_44, ID_45, ID_46, 
                   ID_47, ID_48, ID_49, ID_50, ID_51, ID_52, ID_53, ID_54, ID_55, ID_56, ID_57, ID_58, ID_59, ID_60, ID_61, 
                   ID_62, ID_63, ID_64, ID_65, ID_66, ID_67, ID_68, ID_69, ID_70, ID_71, ID_72, ID_73, ID_74, ID_75, ID_76,
                   ID_77, ID_78, ID_79, ID_80, ID_81, ID_82, ID_83, ID_84, ID_85, ID_86, ID_87, ID_88, ID_89, ID_90, ID_91, 
                   ID_92, ID_93, ID_94, ID_95, ID_96, ID_97, ID_98, ID_99, ID_100, ID_101, ID_102, ID_103, ID_104, ID_105)
colnames(countTable) <- c(list.of.samples[1:105])

list.three = c(rep(c(rep("Unmolested",4), rep("Clean.prick",each=3), rep("M.luteus",each=3), rep("E.coli",each=3), 
           rep("S.mar.type",each=3), rep("E.faecalis",each=3), rep("P.rettgeri",each=3), rep("Ecc15",each=3), 
           "S.aureus", "P.sneebia", "S.mar.DB11","P.entomophila", rep("E.faecalis.heat",each=3), 
           rep("P.rettgeri.heat",each=3)),3))
time.three = c(
    rep(
        c(c("zero"),
        rep(c("twelve","thirty.six","five.half"), 8),
        rep("twelve", 4),
        rep(c("twelve","thirty.six","five.half"),2))
        ,3))

design <- data.frame(row.names = colnames(countTable), 
                     treatment = list.three, time = time.three,
                     rep = c(rep(1:3, each = 35))
                     )

CountTable <- countTable[(1:(dim(countTable)[1]-5)),]

write.table(CountTable, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_count_of_all_samples.txt") #17558 genes
write.table(design, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_experimental_design_for_all_samples.txt")

#Filter the hits with less than 5 reads
filteredCountTable <- vector()
length.of.table <- as.numeric(dim(CountTable)[1])
width.of.table <- as.numeric(dim(CountTable)[2])
indicator <- NULL
for (i in 1:length.of.table){
    for (j in 1:width.of.table){
        #cat("ith row: ", i, "jth column: ",j,"value: ",CountTable[i,j], "\n")
        indicator <- isTRUE(CountTable[i,j] < 5)
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break
    }
    if (j == width.of.table){
        filteredCountTable <- rbind(filteredCountTable,CountTable[i,1:width.of.table])   
    }
}
CountTable <- filteredCountTable
write.table(CountTable, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_count_of_all_samples_more_than_5_reads.txt") #7064 genes

#=====
library("DESeq")

CountTable <- read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_count_of_all_samples_more_than_5_reads.txt")
design <- read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_experimental_design_for_all_samples.txt")

#Sort out CountTable by hour // excluding unmolested and clean prick controls
CountTable.12hr <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,
                                 78,81,84,87,90,93,96,97,98,99,100,103)]
design.12hr <- design[c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,
                    78,81,84,87,90,93,96,97,98,99,100,103),]
CountTable.36hr <- CountTable[,c(9,12,15,18,21,24,31,34,44,47,50,53,56,59,66,69,79,82,85,88,91,94,101,104)]
design.36hr <- design[c(9,12,15,18,21,24,31,34,44,47,50,53,56,59,66,69,79,82,85,88,91,94,101,104),]
CountTable.5.5d <- CountTable[,c(10,13,16,19,22,25,32,35,45,48,51,54,57,60,67,70,80,83,86,89,92,95,102,105)]
design.5.5d <- design[c(10,13,16,19,22,25,32,35,45,48,51,54,57,60,67,70,80,83,86,89,92,95,102,105),]

#Sort out CountTable by treatment // excluding unmolested and clean prick controls
#CountTable.treatment <- CountTable[,c(9:35,43:70,78:105)] #excluding unmolested and clean prick controls
#design.treatment <- design[c(9:35,43:70,78:105),] #These do not work, as "Parametric dispersion fit failed."
CountTable.M.luteus <- CountTable[,c(8,9,10,43,44,45,78,79,80)] #These work only when method="pooled-CR"
design.M.luteus <- design[c(8,9,10,43,44,45,78,79,80),]
CountTable.E.faecalis.with.heat <- CountTable[,c(17,18,19,30,31,32,52,53,54,65,66,67,87,88,89,100,101,102)] #These work only when method="pooled-CR"
design.E.faecalis.with.heat <- design[c(17,18,19,30,31,32,52,53,54,65,66,67,87,88,89,100,101,102),]
CountTable.Ecc15 <- CountTable[,c(23,24,25,58,59,60,93,94,95)] #These work only when method="pooled-CR"
design.Ecc15 <- design[c(23,24,25,58,59,60,93,94,95),]
CountTable.E.coli <- CountTable[,c(11,12,13,46,47,48,81,82,83)] #These work only when method="pooled-CR"
design.E.coli <- design[c(11,12,13,46,47,48,81,82,83),]
CountTable.S.mar.type.and.db11 <- CountTable[,c(14,15,16,28,49,50,51,63,84,85,86,98)] #These work only when method="pooled-CR"
design.S.mar.type.and.db11 <- design[c(14,15,16,28,49,50,51,63,84,85,86,98),]
CountTable.P.rettgeri.with.heat <- CountTable[,c(20,21,22,33,34,35,55,56,57,68,69,70,90,91,92,103,104,105)] #These work only when method="pooled-CR"
design.P.rettgeri.with.heat <- design[c(20,21,22,33,34,35,55,56,57,68,69,70,90,91,92,103,104,105),]
CountTable.P.rettgeri.and.P.sneebia <- CountTable[,c(20,21,22,27,55,56,57,62,90,91,92,97)] #These work only when method="pooled-CR"
design.P.rettgeri.and.P.sneebia <- design[c(20,21,22,27,55,56,57,62,90,91,92,97),]

#Enter your choice // #CountTable.12hr, design.12hr, etc
CountTable <-CountTable.P.rettgeri.and.P.sneebia
design <- design.P.rettgeri.and.P.sneebia

#Normalization: estimate the effective library size.
cds = newCountDataSet(CountTable, design)
cds = estimateSizeFactors(cds)
sizeFactors(cds) #now divide each column of the count table by this size factor
head(counts(cds, normalized=TRUE))

#Diagnostics - PCA
#All the time points
cdsB = estimateDispersions(cds, method ="blind")
#or (ex. M.luteus, Ecc15, E.coli)
cdsB = estimateDispersions(cds, method ="pooled-CR", modelFormula = count ~ time) 
#or (ex. E.faecalis.with.heat, S.mar.type.and.db11, P.rettgeri.with.heat, P.rettgeri.and.P.sneebia)
cdsB = estimateDispersions(cds, method ="pooled-CR", modelFormula = count ~ time + treatment) 
vsd = varianceStabilizingTransformation(cdsB)

#General case: Throwing everything into one PCA
CountTable = CountTable[,c(5:35,40:70,75:105)] #Exclude the unmolested controls = 93 samples
design = design[c(5:35,40:70,75:105),]

cds = newCountDataSet(CountTable, design) #design 
cds = estimateSizeFactors(cds) #Normalization: estimate the effective library size.
sizeFactors(cds)
#cdsB = estimateDispersions(cds, method ="pooled-CR", modelFormula = count ~ time + treatment) #doesn't work "parametic dispersion fit failed" 
#cdsB = estimateDispersions(cds, method ="pooled", modelFormula = count ~ time + treatment) #doesn't work "None of your conditions is replicated"
#cdsB = estimateDispersions(cds, method ="per-condition", modelFormula = count ~ time + treatment) #doesn't work becasue it is a multivariate design
#cdsB = estimateDispersions(cds, method ="blind")

#for DESeq2, 
#cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design, design = ~time + treatment)
#cds = estimateSizeFactors(cds)
#cdsB = estimateDispersions(cds) #for DESeq2 -- this one worked even though it gave some errors in the beginning

vsd = varianceStabilizingTransformation(cdsB)
#plotPCA(vsd, intgroup=c("time","treatment"))
##OR Use this custom code from https://gist.github.com/igordot/8342684:
plotPCAWithSampleNames = function(x, intgroup=c("time","treatment"), ntop=500){
    library(ggplot2)
    library(genefilter)
    library(lattice)
    
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
    fac.time = factor(apply(as.data.frame(colData(x)[, intgroup[1], drop=FALSE]), 1, paste, collapse=" : "))
    print(fac.time)
    fac.treatment = factor(apply(as.data.frame(colData(x)[, intgroup[2], drop=FALSE]), 1, paste, collapse=" : "))

    # colors
    if( nlevels(fac.treatment) >= 10 ) #number of treatment
        colors = c("darkgoldenrod1", "royalblue1","lightblue1","firebrick2","coral1","pink1","darkslateblue",
                   "deepskyblue3","deepskyblue2", "dodgerblue3","coral4","dodgerblue4","lightskyblue")
    
    pca.data = as.data.frame(pca$x)
    #print(head(pca.data))
    # plot
#     xyplot(
#         PC2 ~ PC1, groups=c(fac.treatment, fac.time), data=pca.data, cex=1.5,
#         aspect = "fill",
#         col = colors,
#         pch = levels(fac.time), ###This is all wrong!! (10/14/15)
#         xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
#         ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
#         panel = function(x, y, ...) {
#             panel.xyplot(x, y, ...);
#             ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
#         },
#         main = draw.key(
#             key = list(
#                 columns=2, rect = list(col = colors), text = list(c(levels(fac.treatment), levels(fac.time))),
#                 #rect = list(c(col = colors), c(pch =shape)), 
#                 rep = FALSE)
#             )
#         )

    X.ggplot =  pca.data[,1]; Y.ggplot =  pca.data[,2]; data = as.matrix(cbind(X.ggplot,Y.ggplot))
    rownames(data) = rownames(pca.data)
    print(head(data)); cat("Variance 1: ", variance[1], "Variance 2: ", variance[2])
    write.table(data, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PC1_and_PC2.txt", quote=F, col.names=T,row.names = T)

    
#     qplot(data[,1], data[,2], data = data, color=levels(fac.treatment), shape=levels(fac.time) ) + 
#         geom_point(size=3) + xlab(paste("PC1 (", variance[1], "%)", sep="")) + ylab(paste("PC2 (", variance[2], "%)", sep=""))
#         #geom_text(hjust=-0.2,vjust=-0.2)
}
#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_everything_using_custom_code.pdf", height=10, width=10) 
#plotPCAWithSampleNames(vsd, intgroup=c("time","treatment"), ntop=500)
#dev.off()

#Plotting a PCA alone
PCs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/mega_RNA-seq_DESeq_PC1_and_PC2.txt", header= T) #A table of PC1 and PC2
design = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_experimental_design_for_all_samples.txt") 
design = design[c(5:35,40:70,75:105),]
PCs.with.design = cbind(PCs, design)
labels=rownames(PCs.with.design)
library("ggplot2")
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=5, aes(color = factor(treatment), shape = factor(time))) +
    scale_color_manual(values = c("darkgoldenrod1", "lightblue1","firebrick2","coral1","lightskyblue","pink1","darkorchid4",
    "blue1","deepskyblue3", "royalblue3","coral4","midnightblue","deepskyblue1")) + xlab("PC1 (45.3%)") + ylab("PC2 (13.1%)") + 
    theme(axis.title = element_text(face = "bold", size = 20),legend.title = element_text(size = 20), legend.text = element_text(size=16)) + ylim(-30,20)
#E.coli and Ecc15 are swapped in the graph compared to the legend. So the colors are assigned accordingly.


#_________
#Specific cases:
#Treatment per Time
#12 hours by treatment
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_12_hrs_by_treatment.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("treatment"))
dev.off()
#36 hours by treatment
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_36_hrs_by_treatment.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("treatment"))
dev.off()
#5.5 day by treatment
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_5.5_days_by_treatment.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("treatment"))
dev.off()

#Time per Treatment
#M.luteus
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_M.luteus_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
#E.faecalis vs E.faecalis heat-killed
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_E.faecalis_with_heatkilled_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time", "treatment"))
dev.off()
#Ecc15
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_Ecc15_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
#E.coli
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_E.coli_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
#S.marcescens Type vs S.marcescens Db11
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_S.mar.type_with_Db11_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time", "treatment"))
dev.off()
#P.rettgeri vs P.rettgeri heat-killed
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_P.rettgeri_with_heatkilled_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time", "treatment"))
dev.off()
#P.rettgeri vs P.sneebia
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_P.rettgeri_and_P.sneebia_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time", "treatment"))
dev.off()

#Add gene names to the big table
CountTable <-read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_count_of_all_samples.csv")
list.of.unknown.genes <- as.data.frame(row.names(CountTable))
#list.of.unknown.genes <- read.table(,header=F) #or any file that contains a list of unknown genes by ID (here, I am using the)
colnames(list.of.unknown.genes) <- c("gene_id")

#Re-read full list
full.list <- read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)

#Create a matrix output
output.mx <- matrix(NA, ncol=1, nrow=dim(list.of.unknown.genes)[1])
colnames(output.mx) <-c("gene_name")

for (i in 1:dim(list.of.unknown.genes)[1]){
    gene.to.compare = as.character(list.of.unknown.genes[i,1])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    #checking.against.full.list = as.matrix(checking.against.full.list)
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}

updated.list <- cbind(list.of.unknown.genes, output.mx, CountTable)
write.table(updated.list, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_count_of_all_samples_with_gene_names.csv", quote =F,row.names=F,col.names=T)

#Check agreement among biological replicates (doing this on my personal computer)
count <- read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_more_than_5_reads.csv")
#count <- log2(count) #has NOT been corrected for the differences in the size of libraries (some have more reads than others)
#That would be the reason why the R.squared value is not as good.
total.sum <- sum(count)/105 #This HAS BEEN attempted to correct for the differences in the size of libraries
list.of.multiples <-list()
for (m in 1:105){
    list.of.multiples[m] <- sum(count[,m])/total.sum
    count[m,] <- count[m,]*list.of.multiples[m]
}
count <- log2(count)

library(ggplot2)
#library(gridExtra)
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/agreement.of.bio.reps.pdf") #Printing out graphs does not seem to work (9/14/15)
r.squared.value.table <- matrix(NA, nrow=35, ncol=3, dimnames=list(c(1:35),c("rep1.vs.rep2","rep1.vs.rep3","rep2.vs.rep3")))
for (i in 1:35){ 
    rep1 <- count[,i]; rep2 <- count[,i+35]; rep3 <- count[,i+70]
    #comparison between rep 1 and rep 2 #gotta say 'print(qplot(...))' for ggplot pdf printing
    FC.lm = lm(rep1 ~ rep2)
    r.squared.value.table[i,1] <- summary(FC.lm)$r.squared
    #comparison between rep 1 and rep 3
    FC.lm = lm(rep1 ~ rep3)
    r.squared.value.table[i,2] <- summary(FC.lm)$r.squared
    #comparison between rep 2 and rep 3
    FC.lm = lm(rep2 ~ rep3)
    r.squared.value.table[i,3] <- summary(FC.lm)$r.squared
    
    #Plot the graphs
    #grid.arrange(
    print(qplot(rep1,rep2,asp=1)+ xlab("Rep 1") + ylab("Rep 2") + ggtitle(paste("Sample_",i,"_replicate agreement")) + theme_classic(base_size = 20) + coord_equal(ratio=1) + geom_abline(intercept=0, slope=1, color="red")) #,
    print(qplot(rep1,rep3,asp=1)+ xlab("Rep 1") + ylab("Rep 3") + ggtitle(paste("Sample_",i,"_replicate agreement")) + theme_classic(base_size = 20) + coord_equal(ratio=1) + geom_abline(intercept=0, slope=1, color="red")) #,
    print(qplot(rep2,rep3,asp=1)+ xlab("Rep 2") + ylab("Rep 3") + ggtitle(paste("Sample_",i,"_replicate agreement")) + theme_classic(base_size = 20) + coord_equal(ratio=1) + geom_abline(intercept=0, slope=1, color="red"))
    #)
    
}
dev.off()
write.table(round(r.squared.value.table,digits=3), file ="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/r.squared.value.table.for.mega.RNA-seq.reps.agreement.txt", quote =F,row.names=T, col.names=T)


#When drawing a PCA with EVERYTHING, do unchallenged samples separate out from the rest? Answer: Yes, sort of (12.01.2015)
#12/01/2015 -- on CBSU
rm(list=ls(all=TRUE)) #delete any previous entry
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)
library("DESeq2")
cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design, design = ~time + treatment)
cds = estimateSizeFactors(cds)
cdsB = estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cdsB)
library(RColorBrewer)
library(genefilter)
library(lattice)
plotPCAWithSampleNames = function(x, intgroup=c("time","treatment"), ntop=500) {
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
    
#     # colors
#     if( nlevels(fac) >= 10 )
#         colors = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"), brewer.pal(12,"Set3"), brewer.pal(7,"Dark2"))
#     
#     # plot
#     xyplot(
#         PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
#         aspect = "fill",
#         col = colors,
#         xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
#         ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
#         panel = function(x, y, ...) {
#             panel.xyplot(x, y, ...);
#             ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
#         },
#         main = draw.key(
#             key = list(
#                 rect = list(col = colors),
#                 text = list(levels(fac)),
#                 rep = FALSE
#             )
#         )
#     )
    
    pca.data = as.data.frame(pca$x)
    X.ggplot =  pca.data[,1]; Y.ggplot =  pca.data[,2]; data = as.matrix(cbind(X.ggplot,Y.ggplot))
    rownames(data) = rownames(pca.data)
    print(head(data)); cat("Variance 1: ", variance[1], "Variance 2: ", variance[2])
    write.table(data, file="mega_RNA-seq_DESeq2_PC1_and_PC2_recalculated_Apr2017.txt", quote=F, col.names=T,row.names = T)
    #write.table(data, file="/local/workdir/ji72/mega_RNA-seq_DESeq2_PC1_and_PC2.txt", quote=F, col.names=T,row.names = T)
}
#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_all_PCA_Dec_2015_predefined.pdf", height=10, width=10) 
#plotPCAWithSampleNames(vsd, intgroup=c("time","treatment"), ntop=500)
#dev.off()

PCs = read.table("/local/workdir/ji72/mega_RNA-seq_DESeq2_PC1_and_PC2.txt", header= T) #A table of PC1 and PC2
PCs.with.design = cbind(PCs, design)
labels=rownames(PCs.with.design)
library("ggplot2")
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_all_PCA_Dec_2015.pdf", height=10, width=10) 
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=4, aes(color = factor(treatment), shape = factor(time))) +
    scale_color_manual(values = c("darkgoldenrod1", "lightblue1","firebrick2","coral1","lightskyblue","pink1","darkorchid4","blue1","deepskyblue3", "royalblue3","coral4","midnightblue","deepskyblue1","forestgreen")) + xlab("PC1 (46.1%)") + ylab("PC2 (12.7%)") + 
    theme(axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=16)) + ylim(-30,20)
ggplot.PC
dev.off()


##Drawing a PCA with 12hr samples
#12/17/2015 -- on CBSU
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)

library("DESeq2")
#Only choosing 12hr samples (both live and heatkilled, but not clean prick and unchallenged)
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103),]
# design <- design[c(8,11,14,17,20,23,26,27,28,29,43,46,49,52,55,58,61,62,63,64,78,81,84,87,90,93,96,97,98,99),]
attach(design)

cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design, design = ~ treatment)
cds = estimateSizeFactors(cds)
cdsB = estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cdsB)

pdf(file="Rplot.pdf")
plotPCA(vsd, intgroup=c("treatment", "rep"))
dev.off()

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
    print(head(data)); cat("Variance 1: ", variance[1], "Variance 2: ", variance[2])
    write.table(data, file="/local/workdir/ji72/mega_RNA-seq_DESeq2_12hr_only_PC1_and_PC2.txt", quote=F, col.names=T,row.names = T)
}
plotPCAWithSampleNames(vsd, intgroup=c("treatment"), ntop=500)

#PCs = read.table("/local/workdir/ji72/mega_RNA-seq_DESeq2_12hr_only_PC1_and_PC2.txt", header= T) #A table of PC1 and PC2

#Change the treatment names to something more audience-friendly
treat = as.character(design$treatment)
treat=sub("S.mar.DB11", "S.marcescens Db11",treat); treat=sub("S.mar.type", "S.marcescens Type",treat)
treat=sub("E.faecalis.heat", "E.faecalis heatkilled",treat); treat=sub("P.rettgeri.heat", "P.rettgeri heatkilled",treat)
treat=sub("P. entomophila", "P.entomophila",treat)
design = cbind(treat,design[,2:3]); colnames(design) = c("treatment","time","rep")
attach(design)
PCs.with.design = cbind(PCs, design)
labels=rownames(PCs.with.design)

#modified on 4/7/2016
library("ggplot2")
cols = c("M.luteus" = "lightpink", "E.coli" = "darkseagreen2", "S.marcescens Type" = "lightseagreen", "Ecc15" = "lightskyblue", 
         "P.rettgeri" = "dodgerblue2","E.faecalis" = "coral1", "S.aureus" = "orangered3","P.sneebia" = "blue",
         "S.marcescens Db11" = "navy", "P.entomophila" = "purple3", "E.faecalis heatkilled" = "dark salmon","P.rettgeri heatkilled" = "steelblue1")

#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_12hr_only_PCA_Dec_2015.pdf", height=10, width=13) #used different colors
#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_12hr_only_PCA_Apr_2016.pdf", height=10, width=13)
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq2_12hr_live_only_PCA_Apr_2016.pdf", height=10, width=13)
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=8, aes(color = factor(treatment))) +
    scale_color_manual(values = cols, breaks=c("M.luteus", "E.coli", "S.marcescens Type","Ecc15","P.rettgeri","P.rettgeri heatkilled","E.faecalis","E.faecalis heatkilled","S.aureus","P.sneebia","S.marcescens Db11","P.entomophila"), name="Infection Conditions") 
+ xlab("PC1 (39.3%)") + ylab("PC2 (23.2%)") + 
#+ xlab("PC1 (34.0%)") + ylab("PC2 (27.3%)") + 
    theme_bw() + theme(axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=20), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) #+ xlim(-22,20) + ylim(-22,20)

ggplot.PC
dev.off()

#4/7/2016, 12hr live + heatkilled, PC1: 34%, PC2: 27.3%
#4/20/2016, 12hr live, PC1: 39.3%, PC2: 23.2% 

col2rgb("darkorchid4", alpha = FALSE)


#2/1: change colors of the dots 
ggplot.PC = ggplot(PCs.with.design, aes(x=X.ggplot, y=Y.ggplot)) + geom_point(size=7, aes(color = factor(treatment))) + scale_color_manual(values=cols, breaks=c("Uninjured", "Sterile wounded",
                                                                                                                                                                 "M.luteus", "E.coli", "S.marcescens Type","Ecc15","P.rettgeri","E.faecalis","S.aureus","P.sneebia","S.marcescens Db11","P.entomophila"), name="Infection Condition") + xlab("PC1 (40.0%)") + ylab("PC2 (21.9%)") +
    theme_bw() + theme(legend.justification=c(0,0), legend.position=c(0,0), axis.title = element_text(face = "bold", size = 20), legend.title = element_text(size = 20), legend.text = element_text(size=20), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) #+ xlim(-20,20) + ylim(-20,20)
ggplot.PC

