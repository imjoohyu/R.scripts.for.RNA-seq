#DESeq Commands for Mega RNA-seq for basic comparison
#Date: September 16, 2015
#ji72

rm(list=ls(all=TRUE))

library("DESeq")
#Read in the mega table
setwd("/workdir/ji72")
#setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
#CountTable <- read.table("mega_RNA-seq_count_of_all_samples_more_than_5_reads.txt")
#design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt")
CountTable <- read.table("mega_RNA-seq_count_of_all_samples_more_than_5_reads.txt")
#design <- read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_experimental_design_for_all_samples.txt")

#Breakdown of the samples
#Unchallenged.averaged
Unchallenged.averaged = matrix(NA, nrow=dim(CountTable)[1], ncol=3)
for (i in 1:dim(CountTable)[1]){
    Unchallenged.averaged[i,1] <- sum(CountTable[i,2],CountTable[i,3],CountTable[i,4])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,2] <- sum(CountTable[i,37],CountTable[i,38],CountTable[i,39])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,3] <- sum(CountTable[i,72],CountTable[i,73],CountTable[i,74])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
}
colnames(Unchallenged.averaged) <-c("UC.Rep1", "UC.Rep2", "UC.Rep3")
rownames(Unchallenged.averaged) <- row.names(CountTable)
Unchallenged.averaged <- round(Unchallenged.averaged, digits=0)

#clean.prick
clean.prick.12hr = CountTable[,c(5,40,75)]; clean.prick.36hr = CountTable[,c(6,41,76)]; clean.prick.5.5d = CountTable[,c(7,42,77)]
#M.luteus
M.luteus.12hr = CountTable[,c(8,43,78)]; M.luteus.36hr = CountTable[,c(9,44,79)]; M.luteus.5.5d = CountTable[,c(10,45,80)]
#E.coli
E.coli.12hr = CountTable[,c(11,46,81)] ; E.coli.36hr = CountTable[,c(12,47,82)]; E.coli.5.5d = CountTable[,c(13,48,83)]
#S.mar.type
S.mar.type.12hr = CountTable[,c(14,49,84)]; S.mar.type.36hr = CountTable[,c(15,50,85)] ; S.mar.type.5.5d = CountTable[,c(16,51,86)]
#E.fae.live
E.fae.live.12hr = CountTable[,c(17,52,87)]; E.fae.live.36hr = CountTable[,c(18,53,88)]; E.fae.live.5.5d = CountTable[,c(19,54,89)];
#P.rett.live
P.rett.live.12hr = CountTable[,c(20,55,90)]; P.rett.live.36hr = CountTable[,c(21,56,91)]; P.rett.live.5.5d = CountTable[,c(22,57,92)];
#Ecc15
Ecc15.12hr = CountTable[,c(23,58,93)]; Ecc15.36hr = CountTable[,c(24,59,94)]; Ecc15.5.5d = CountTable[,c(25,60,95)]
#Virulent ones
S.aureus.12hr = CountTable[,c(26,61,96)] ; P.sneebia.12hr = CountTable[,c(27,62,97)]; S.mar.Db11.12hr = CountTable[,c(28,63,98)]; P.ento.12hr = CountTable[,c(29,64,99)]
#E.fae.heatkilled
E.fae.heatkilled.12hr = CountTable[,c(30,65,100)]; E.fae.heatkilled.36hr = CountTable[,c(31,66,101)]; E.fae.heatkilled.5.5d = CountTable[,c(32,67,102)];
#P.rett.heat
P.rett.heatkilled.12hr = CountTable[,c(33,68,103)]; P.rett.heatkilled.36hr = CountTable[,c(34,69,104)]; P.rett.heatkilled.5.5d = CountTable[,c(35,70,105)];

#list of conditions
list.of.conditions <- list(clean.prick.12hr, clean.prick.36hr, clean.prick.5.5d, M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d,
                           S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d, P.rett.live.12hr, P.rett.live.36hr,
                           P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr,
                           E.fae.heatkilled.12hr, E.fae.heatkilled.36hr, E.fae.heatkilled.5.5d, P.rett.heatkilled.12hr, P.rett.heatkilled.36hr, P.rett.heatkilled.5.5d)
list.of.name.of.conditions <- list("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d",
                           "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                           "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
                           "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")

res.of.all.conditions <- matrix(rownames(CountTable), nrow=dim(CountTable)[1])

for (j in 1:length(list.of.conditions)){
    
    cat("We are working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
    condition.count.table <- as.matrix(data.frame(list.of.conditions[j]))
    count.table <- cbind(Unchallenged.averaged, condition.count.table)
    design.sp = data.frame(row.names = colnames(count.table), condition = c("Unchallenged","Unchallenged","Unchallenged",
                as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]),as.character(list.of.name.of.conditions[j])) 
                )#,rep = c("1","2","3","1","2","3"))

    cds = newCountDataSet(count.table, design.sp)
    #Normalization: estimate the effective library size.
    cds = estimateSizeFactors(cds)
    sizeFactors(cds) #now divide each column of the count table by this size factor

    #Estimate variance using the CR-adjusted profile likelihood dispersion estimates
    cds = estimateDispersions(cds, method = "pooled-CR", modelFormula = count ~ condition)
    #Inspect the estimated dispersions
    #pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_estimate_of_dispersion_pooled.pdf", height=5, width=5)

    res = nbinomTest(cds, "Unchallenged", as.character(list.of.name.of.conditions[j])) #the last twos are condition names.
    
    #Plot a dispersion plot, an MA plot and a p-value plot for a given condition
    png(file = paste("Mega_RNA-seq_dispersion_MA_plot_and_pval_histogram_of_",list.of.name.of.conditions[j],".jpg", sep=""), width=1300, height=400)
    par(mfrow = c(1,3))
    plotDispEsts(cds, cex=1)
    plotMA(res, cex=1) 
    hist(res$pval, breaks=100, col="skyblue", freq=F, border="slateblue", 
        main=paste("hist(p-val) between 
                   (Averaged) unchallenged vs", 
                   list.of.name.of.conditions[j]), sep="")
    dev.off()

    #Save the output
    res.with.specific.info <- cbind(res[,6], res[,8]); colnames(res.with.specific.info) <- c("log2FC","padj") #log2 Fold change and padj
    res.of.all.conditions <- cbind(res.of.all.conditions, res.with.specific.info)
    
    cat("Finished! We are done working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
}

write.table(res.of.all.conditions, file="Mega_RNA-seq_DESeq_basic_comparison_of_conditions_to_averaged_unchallenged_fold_change.txt", quote=F, row.names=F, col.names=T)





#Sig = res[which(res$padj < 0.1),] #cutoff: FDR of 10%
#resSig <- resSig[order(resSig$log2FoldChange, decreasing=T),]
#ten.most.upregulated.0.8r <- head(resSig[order(resSig$log2FoldChange, decreasing=T),],10) #10 most upregulated
#ten.most.downregulated.0.8r <- head(resSig[order(resSig$log2FoldChange, decreasing=F),],10)#10 most downregulated
#number.of.genes.with.significant.differential.expression.0.8r <- table(res$padj < 0.1) # shows the number of genes with significant differential expression at a FDR of 10%
