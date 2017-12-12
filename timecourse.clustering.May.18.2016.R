#Option 3: Only use the genes that were DE in all 3 time points (12hr, 36hr, and 5.5d samples).
#DE between each of these time points and the 0hr unchallenged samples. This way may show a stronger impact of infection
#After pulling out the DE genes, draw the clustering heatmap using the [EBSeq-HMM] method.
rm(list=ls(all=TRUE))
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)

#(i) Separate the table into 10 different conditions and pull out only the genes that are DE in all time points
de.list = de.list[,c(1,2,9:52)]
M.luteus = de.list[,c(1:2,3:8)]; E.coli = de.list[,c(1:2,9:14)]; S.mar.type = de.list[,c(1:2,15:20)]; E.fae.live = de.list[,c(1:2,21:26)]; P.rett.live =de.list[,c(1:2,27:32)]; Ecc15 = de.list[,c(1:2,33:38)]

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

for (i in 1:length(list.of.bacteria.name)){
    condition.count.table = list.of.bacteria[[i]]
    sig.DE.only = condition.count.table[!grepl("N",condition.count.table[,4]),]
    sig.DE.only = sig.DE.only[!grepl("N",sig.DE.only[,6]),]
    sig.DE.only = sig.DE.only[!grepl("N",sig.DE.only[,8]),]
    write.table(sig.DE.only, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[i],"_sig_DE_in_all_time_points.txt",sep=""), quote=F, row.names = F, col.names = T)
}

#(ii) Here's the list of DE genes that are significant in all time points in each infection condition
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_all_time_points.txt", header=T)
M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_all_time_points.txt", header=T)
E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_all_time_points.txt", header=T)
S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_all_time_points.txt", header=T)
E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_all_time_points.txt", header=T)
P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_all_time_points.txt", header=T)
Ecc15.gene.only = Ecc15.gene.only[,1]

#(iii) Perform EBSeq-HMM on the local machine
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
CountTable = CountTable[,c(2:35, 37:70, 72:105)] #excluding unchallenged 0hr samples, 17558 x 102
colnames(CountTable) = c("UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
                         "M.luteus.12hr.Rep1","M.luteus.36hr.Rep1","M.luteus.5.5d.Rep1","E.coli.12hr.Rep1", "E.coli.36hr.Rep1","E.coli.5.5d.Rep1",
                         "S.mar.type.12hr.Rep1","S.mar.type.36hr.Rep1","S.mar.type.5.5d.Rep1","E.fae.live.12hr.Rep1","E.fae.live.36hr.Rep1","E.fae.live.5.5d.Rep1",
                         "P.rett.live.12hr.Rep1","P.rett.live.36hr.Rep1","P.rett.live.5.5d.Rep1","Ecc15.12hr.Rep1","Ecc15.36hr.Rep1","Ecc15.5.5d.Rep1",
                         "S.aureus.12hr.Rep1","P.sneebia.12hr.Rep1","S.mar.Db11.12hr.Rep1","P.ento.12hr.Rep1","E.fae.heatkilled.12hr.Rep1","E.fae.heatkilled.36hr.Rep1","E.fae.heatkilled.5.5d.Rep1",
                         "P.rett.heatkilled.12hr.Rep1","P.rett.heatkilled.36hr.Rep1","P.rett.heatkilled.5.5d.Rep1","UC.Rep4","UC.Rep5","UC.Rep6",
                         "clean.prick.12hr.Rep2","clean.prick.36hr.Rep2","clean.prick.5.5d.Rep2",
                         "M.luteus.12hr.Rep2","M.luteus.36hr.Rep2","M.luteus.5.5d.Rep2","E.coli.12hr.Rep2", "E.coli.36hr.Rep2","E.coli.5.5d.Rep2",
                         "S.mar.type.12hr.Rep2","S.mar.type.36hr.Rep2","S.mar.type.5.5d.Rep2","E.fae.live.12hr.Rep2","E.fae.live.36hr.Rep2","E.fae.live.5.5d.Rep2",
                         "P.rett.live.12hr.Rep2","P.rett.live.36hr.Rep2","P.rett.live.5.5d.Rep2","Ecc15.12hr.Rep2","Ecc15.36hr.Rep2","Ecc15.5.5d.Rep2",
                         "S.aureus.12hr.Rep2","P.sneebia.12hr.Rep2","S.mar.Db11.12hr.Rep2","P.ento.12hr.Rep2","E.fae.heatkilled.12hr.Rep2","E.fae.heatkilled.36hr.Rep2","E.fae.heatkilled.5.5d.Rep2",
                         "P.rett.heatkilled.12hr.Rep2","P.rett.heatkilled.36hr.Rep2","P.rett.heatkilled.5.5d.Rep2",
                         "UC.Rep7","UC.Rep8","UC.Rep9","clean.prick.12hr.Rep3","clean.prick.36hr.Rep3","clean.prick.5.5d.Rep3",
                         "M.luteus.12hr.Rep3","M.luteus.36hr.Rep3","M.luteus.5.5d.Rep3","E.coli.12hr.Rep3", "E.coli.36hr.Rep3","E.coli.5.5d.Rep3",
                         "S.mar.type.12hr.Rep3","S.mar.type.36hr.Rep3","S.mar.type.5.5d.Rep3","E.fae.live.12hr.Rep3","E.fae.live.36hr.Rep3","E.fae.live.5.5d.Rep3",
                         "P.rett.live.12hr.Rep3","P.rett.live.36hr.Rep3","P.rett.live.5.5d.Rep3","Ecc15.12hr.Rep3","Ecc15.36hr.Rep3","Ecc15.5.5d.Rep3",
                         "S.aureus.12hr.Rep3","P.sneebia.12hr.Rep3","S.mar.Db11.12hr.Rep3","P.ento.12hr.Rep3","E.fae.heatkilled.12hr.Rep3","E.fae.heatkilled.36hr.Rep3","E.fae.heatkilled.5.5d.Rep3",
                         "P.rett.heatkilled.12hr.Rep3","P.rett.heatkilled.36hr.Rep3","P.rett.heatkilled.5.5d.Rep3")

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d); M.luteus = M.luteus[row.names(M.luteus) %in% M.luteus.gene.only,]
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d); E.coli = E.coli[row.names(E.coli) %in% E.coli.gene.only,]
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d); S.mar.type = S.mar.type[row.names(S.mar.type) %in% S.mar.type.gene.only,]
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d); E.fae.live = E.fae.live[row.names(E.fae.live) %in% E.fae.live.gene.only,]
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d); P.rett.live = P.rett.live[row.names(P.rett.live) %in% P.rett.live.gene.only,]
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d); Ecc15 = Ecc15[row.names(Ecc15) %in% Ecc15.gene.only,]

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

#Perform EBSeq-HMM
for (j in 1:length(list.of.bacteria.name)){
    cat("EBSeq-HMM - We are working on the following comparison: 0hr-12hr-36hr-5.5d of ", as.character(list.of.bacteria.name[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 2){ #getting rid of Rep2 of E.coli 36hr
        CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d")} #removing 36hr rep2 sample 
    else {CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d") }    
    Conditions <- factor(CondVector, levels=c("t.0hr","t.12hr","t.36hr","t.5.5d"))
    
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=count.table, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
    EBSeqHMMGeneOut$MgAllMAPChar #classifies the pattern of gene expression over time per gene
    AllPaths = GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=TRUE) #Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
    gene.list.for.a.given.path = GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, cutoff=0.8, OnlyDynamic=TRUE,Paths=NULL) #Obtain confident gene calls for classifying genes into expression paths, PP>=0.5 is most likely path
    number.of.genes.in.each.category = gene.list.for.a.given.path$NumEach; number.of.genes.in.each.category #already sorted
    genes.in.each.category = gene.list.for.a.given.path$EachPathNames; genes.in.each.category
    de.list = GetDECalls(EBSeqHMMGeneOut,FDR=.05) #Another way to obtain DE gene/isoform list at a certain FDR (but not PPs)
    
    #For a specific gene, draw a plot of expression
    #PlotExp(count.table, Conditions, "Gene_1")
    
    write.table(de.list, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM/EBSeqHMM_DE_results_for_",list.of.bacteria.name[j],".txt", sep=""), quote=F)
}

setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
count.table = read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)

#Do the case for M.luteus
s=1; most.likely.paths = read.table(paste("clustering/EBSeq-HMM/EBSeqHMM_DE_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)

exp.list = c()
for (n in 1:dim(most.likely.paths)[1]){
    gene.to.search = row.names(most.likely.paths[n,])
    gene.found = M.luteus[which(row.names(M.luteus) == gene.to.search),]
    exp.list = rbind(exp.list, gene.found)
}
mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:9],1,mean),apply(exp.list[,10:12],1,mean))
mean.bacteria.log10 = log10(mean.bacteria + 1)
most.likely.paths.M.luteus = cbind(most.likely.paths, mean.bacteria.log10)
colnames(most.likely.paths.M.luteus) = c("paths","max_pp","0hr","12hr","36hr","5.5d")
M.luteus.reshape = reshape(most.likely.paths.M.luteus, direction="long", varying = 3:6, idvar = "ids", ids = rownames(most.likely.paths.M.luteus), timevar="time", v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))

#12 patterns
pattern = levels(M.luteus.reshape$paths)
for (i in 1:length(levels(M.luteus.reshape$paths))){
    assign(paste0("M.luteus.reshape.c",i), M.luteus.reshape[which(M.luteus.reshape$paths == pattern[i]),])
}
library(ggplot2); library(grid); library(gridExtra)
q1 = ggplot(data = M.luteus.reshape.c1, aes(x=M.luteus.reshape.c1$time, y = M.luteus.reshape.c1$log10.exp.value, group = M.luteus.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = M.luteus.reshape.c2, aes(x=M.luteus.reshape.c2$time, y = M.luteus.reshape.c2$log10.exp.value, group = M.luteus.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = M.luteus.reshape.c3, aes(x=M.luteus.reshape.c3$time, y = M.luteus.reshape.c3$log10.exp.value, group = M.luteus.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = M.luteus.reshape.c4, aes(x=M.luteus.reshape.c4$time, y = M.luteus.reshape.c4$log10.exp.value, group = M.luteus.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = M.luteus.reshape.c5, aes(x=M.luteus.reshape.c5$time, y = M.luteus.reshape.c5$log10.exp.value, group = M.luteus.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = M.luteus.reshape.c6, aes(x=M.luteus.reshape.c6$time, y = M.luteus.reshape.c6$log10.exp.value, group = M.luteus.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = M.luteus.reshape.c7, aes(x=M.luteus.reshape.c7$time, y = M.luteus.reshape.c7$log10.exp.value, group = M.luteus.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = M.luteus.reshape.c8, aes(x=M.luteus.reshape.c8$time, y = M.luteus.reshape.c8$log10.exp.value, group = M.luteus.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = M.luteus.reshape.c9, aes(x=M.luteus.reshape.c9$time, y = M.luteus.reshape.c9$log10.exp.value, group = M.luteus.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = M.luteus.reshape.c10, aes(x=M.luteus.reshape.c10$time, y = M.luteus.reshape.c10$log10.exp.value, group = M.luteus.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = M.luteus.reshape.c11, aes(x=M.luteus.reshape.c11$time, y = M.luteus.reshape.c11$log10.exp.value, group = M.luteus.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = M.luteus.reshape.c12, aes(x=M.luteus.reshape.c12$time, y = M.luteus.reshape.c12$log10.exp.value, group = M.luteus.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q13 = ggplot(data = M.luteus.reshape.c13, aes(x=M.luteus.reshape.c13$time, y = M.luteus.reshape.c13$log10.exp.value, group = M.luteus.reshape.c13$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q14 = ggplot(data = M.luteus.reshape.c14, aes(x=M.luteus.reshape.c14$time, y = M.luteus.reshape.c14$log10.exp.value, group = M.luteus.reshape.c14$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q15 = ggplot(data = M.luteus.reshape.c15, aes(x=M.luteus.reshape.c15$time, y = M.luteus.reshape.c15$log10.exp.value, group = M.luteus.reshape.c15$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15, ncol=5, newpage=T)

grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)
