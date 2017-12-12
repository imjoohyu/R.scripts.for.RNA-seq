#EBSeq-HMM Commands for Mega RNA-seq for basic comparison
#Identifies DEGs and specify gene-specific expression paths. Good for ordered conditions.
#Date: September 29, 2015
#ji72

#Done in my laptop
rm(list=ls(all=TRUE))
library("EBSeqHMM")

#Read in the mega table
#setwd("/workdir/edgeR_mega_RNA-seq_Sept_2015")
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable <- read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt") #from deseq.for.samples.with.time.pca.for.mega.rna.R (a simple sorting >5 reads code that has nothing to do with deseq yielded this result)
#design <- read.table("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_experimental_design_for_all_samples.txt")

#Label the samples accordingly
Unchallenged.averaged = matrix(NA, nrow=dim(CountTable)[1], ncol=3)
for (i in 1:dim(CountTable)[1]){
    Unchallenged.averaged[i,1] <- sum(CountTable[i,2],CountTable[i,3],CountTable[i,4])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,2] <- sum(CountTable[i,37],CountTable[i,38],CountTable[i,39])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,3] <- sum(CountTable[i,72],CountTable[i,73],CountTable[i,74])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
}
colnames(Unchallenged.averaged) <-c("UC.Rep1", "UC.Rep2", "UC.Rep3")
rownames(Unchallenged.averaged) <- row.names(CountTable)
Unchallenged.averaged <- round(Unchallenged.averaged, digits=0)

clean.prick = CountTable[,c(5,40,75,6,41,76,7,42,77)] #clean.prick
M.luteus = CountTable[,c(8,43,78,9,44,79,10,45,80)] #M.luteus
E.coli = CountTable[,c(11,46,81,12,82,13,48,83)] #E.coli
S.mar.type = CountTable[,c(14,49,84,15,50,85,16,51,86)] #S.mar.type
E.fae.live = CountTable[,c(17,52,87,18,53,88,19,54,89)] #E.fae.live
P.rett.live = CountTable[,c(20,55,90,21,56,91,22,57,92)] #P.rett.live
Ecc15 = CountTable[,c(23,58,93,24,59,94,25,60,95)] #Ecc15
E.fae.heatkilled = CountTable[,c(30,65,100,31,66,101,32,67,102)] #E.fae.heatkilled
P.rett.heatkilled = CountTable[,c(33,68,103,34,69,104,35,70,105)] #P.rett.heatkilled
S.aureus.12hr = CountTable[,c(26,61,96)]; P.sneebia.12hr = CountTable[,c(27,62,97)]; S.mar.Db11.12hr = CountTable[,c(28,63,98)]; P.ento.12hr = CountTable[,c(29,64,99)]

#list of conditions
list.of.conditions <- list(clean.prick, M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live,
                           Ecc15, E.fae.heatkilled, P.rett.heatkilled) #S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr
list.of.name.of.conditions <- list("clean.prick", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", 
                                   "E.fae.heatkilled", "P.rett.heatkilled")  #"S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr"
de.of.all.conditions <- matrix(rownames(CountTable), nrow=dim(CountTable)[1], ncol=1)

#Execute EBSeq-HMM on 12hr, 36hr, and 5.5d time points. I excluded 0hr because it messes up the first point.

for (j in 7:7){ #length(list.of.conditions
    cat("EBSeq-HMM - We are working on the following comparison: 12hr-36hr-5.5d of ", as.character(list.of.name.of.conditions[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.conditions[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 3){ #getting rid of Rep2 of E.coli 36hr
        CondVector <- c("t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d")} #removing 36hr rep2 sample 
    else {CondVector <- c("t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d") }    
    Conditions <- factor(CondVector, levels=c("t.12hr","t.36hr","t.5.5d"))
    
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=count.table, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
    EBSeqHMMGeneOut$MgAllMAPChar #classifies the pattern of gene expression over time per gene
    AllPaths = GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=TRUE) #Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
    gene.list.for.a.given.path = GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, cutoff=0.8, OnlyDynamic=TRUE,Paths=NULL) #Obtain confident gene calls for classifying genes into expression paths, PP>=0.5 is most likely path
    number.of.genes.in.each.category = gene.list.for.a.given.path$NumEach; number.of.genes.in.each.category #already sorted
    genes.in.each.category = gene.list.for.a.given.path$EachPathNames; genes.in.each.category
    de.list = GetDECalls(EBSeqHMMGeneOut,FDR=.05) #Another way to obtain DE gene/isoform list at a certain FDR (but not PPs)
    
    #For a specific gene, draw a plot of expression
    #PlotExp(count.table, Conditions, "Gene_1")
    
}
