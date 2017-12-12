#GeoDE Commands for Mega RNA-seq for basic comparison
#GeoDE identifies DEGs in a multivariate way by using Characteristid Directions
#unlike DESeq and edgeR (univariate) so it accounts for gene-gene statistical dependencies.
#Date: September 28, 2015
#ji72

#Done in my laptop
rm(list=ls(all=TRUE))
library("GeoDE")

#Read in the mega table
#setwd("/workdir/edgeR_mega_RNA-seq_Sept_2015")
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable <- read.table("mega_RNA-seq_count_of_all_samples_more_than_5_reads.txt") #from deseq.for.samples.with.time.pca.for.mega.rna.R (a simple sorting >5 reads code that has nothing to do with deseq yielded this result)
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

#clean.prick
clean.prick.12hr = CountTable[,c(5,40,75)]; clean.prick.36hr = CountTable[,c(6,41,76)]; clean.prick.5.5d = CountTable[,c(7,42,77)]
#M.luteus
M.luteus.12hr = CountTable[,c(8,43,78)]; M.luteus.36hr = CountTable[,c(9,44,79)]; M.luteus.5.5d = CountTable[,c(10,45,80)]
#E.coli
E.coli.12hr = CountTable[,c(11,46,81)] ; E.coli.36hr = CountTable[,c(12,82)]; E.coli.5.5d = CountTable[,c(13,48,83)]
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

de.of.all.conditions <- matrix(rownames(CountTable), nrow=dim(CountTable)[1], ncol=1)

#Execute GeoDE
gammas =c(1.0) #shrinkage value
data(AllGMTfiles) #data(GeneOntology_BP.gmt) #Load the GMT file
for (j in 25:25){ #length(list.of.conditions
    cat("GeoDE - We are working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
    condition.count.table <- as.matrix(data.frame(list.of.conditions[j]))
    count.table <- cbind(Unchallenged.averaged, condition.count.table)
    if (j == 8){ #for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr
        design.sp = as.factor(c(1,1,1,2,2))} #removing 36hr rep2 sample 
    else{ design.sp = as.factor(c(1,1,1,2,2,2)) }    
    count.table = cbind(as.data.frame(rownames(count.table)), count.table); row.names(count.table) = NULL
    
    chdir_analysis <- chdirAnalysis(count.table,design.sp,gammas,CalculateSig=TRUE,nnull=10) #Characteristic Direction Analysis
    #lapply(chdir_analysis$chdirprops[[1]],head) #gives the CD vectors
    top.ten = lapply(chdir_analysis$results, function(x) x[1:10]) #show the first few of the most important genes.
    print(top.ten); coef = chdir_analysis$chdirprops$chdir
    num.sig.genes = chdir_analysis$chdirprops$number_sig_genes
    print("The estimated number of significant genes: ")
    print(num.sig.genes)
    
    #Run PAEA test
    #chdir_analysis <- chdirAnalysis(count.table,design.sp,gammas,CalculateSig=F)
    #PAEAtest = multigmtPAEAAnalysis(chdir_analysis$chdirprops, AllGMTfiles[6], gammas) #Run the PAEAtest with GO, GeneSigDB, KEGG
}

#Comments: 9/28/15, 11pm
#Main question: How do I know that they are significant??? and How does shrinkage value work?
