#################################################################################################################
#Use EBSeq-HMM to make clusters of DE genes
#Modified on 5/22/2016
#Joo Hyun Im (ji72)

#Only use the genes that were DE in all 3 time points (12hr, 36hr, and 5.5d samples).
#DE between each of these time points and the 0hr unchallenged samples. This way may show a stronger impact of infection
#After pulling out the DE genes, draw the clustering heatmap using the [EBSeq-HMM] method.

#Done in my laptop
rm(list=ls(all=TRUE))

#(i) Separate the table into 10 different conditions and pull out only the genes that are DE between all time points
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
de.list = de.list[,c(1,2,9:52)]; M.luteus = de.list[,c(1:2,3:8)]; E.coli = de.list[,c(1:2,9:14)]; S.mar.type = de.list[,c(1:2,15:20)]; E.fae.live = de.list[,c(1:2,21:26)]; P.rett.live =de.list[,c(1:2,27:32)]; Ecc15 = de.list[,c(1:2,33:38)]
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
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_all_time_points.txt", header=T); M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_all_time_points.txt", header=T); E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_all_time_points.txt", header=T); S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_all_time_points.txt", header=T); E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_all_time_points.txt", header=T); P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_all_time_points.txt", header=T); Ecc15.gene.only = Ecc15.gene.only[,1]

#(iii) Pull out the expression values of DE genes from the count table
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

#(iv) Perform EBSeq-HMM on the expression of DE genes
library(EBSeqHMM)
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
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

#(v) Draw representation plots
count.table = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
library(ggplot2); library(grid); library(gridExtra)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM/")
for (s in 1:length(list.of.bacteria.name)){ 
    #For a given gene in the most_likely_path, identify the expression values in count.table and cbind
    most.likely.paths = read.table(paste("EBSeqHMM_DE_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
    #most.likely.paths = most.likely.paths[which(most.likely.paths$Max_PP > 0.55),] #getting rid of genes that have not so great categorization
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[s]))
    exp.list = c()
    for (n in 1:dim(most.likely.paths)[1]){
        gene.to.search = row.names(most.likely.paths[n,])
        gene.found = condition.count.table[which(row.names(condition.count.table) == gene.to.search),]
        exp.list = rbind(exp.list, gene.found)
    }
    if (s == 2){
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:8],1,mean),apply(exp.list[,9:11],1,mean))
    }
    else{
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:9],1,mean),apply(exp.list[,10:12],1,mean))
    }
    mean.bacteria.log10 = log10(mean.bacteria + 1) #This could be changed to centered mean
    most.likely.paths.bacteria = cbind(most.likely.paths,mean.bacteria.log10)
    colnames(most.likely.paths.bacteria) = c("paths","max_pp","0hr","12hr","36hr","5.5d")
    most.likely.paths.bacteria.reshape = reshape(most.likely.paths.bacteria, direction="long", varying = 3:6, idvar = "ids", ids = rownames(most.likely.paths.bacteria), timevar="time", v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
    most.likely.paths.bacteria = most.likely.paths.bacteria[order(most.likely.paths.bacteria$paths),]
    #write.table(most.likely.paths.bacteria, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[s],".with.most.liekly.paths.txt",sep=""), quote=F, row.names =T, col.names=T)
    write.table(most.likely.paths.bacteria, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[s],".with.most.liekly.paths.based.on.VTS.txt",sep=""), quote=F, row.names =T, col.names=T)
    
    #Number of patterns
    pattern = levels(most.likely.paths.bacteria.reshape$paths)
    for (i in 1:length(pattern)){
        assign(paste0(list.of.bacteria.name[s],".reshape.c",i), most.likely.paths.bacteria.reshape[which(most.likely.paths.bacteria.reshape$paths == pattern[i]),])
    }
}

#Plot them here
#M.luteus - 12 clusters
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
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#E.coli - 7 clusters
q1 = ggplot(data = E.coli.reshape.c1, aes(x=E.coli.reshape.c1$time, y = E.coli.reshape.c1$log10.exp.value, group = E.coli.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.coli.reshape.c2, aes(x=E.coli.reshape.c2$time, y = E.coli.reshape.c2$log10.exp.value, group = E.coli.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.coli.reshape.c3, aes(x=E.coli.reshape.c3$time, y = E.coli.reshape.c3$log10.exp.value, group = E.coli.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.coli.reshape.c4, aes(x=E.coli.reshape.c4$time, y = E.coli.reshape.c4$log10.exp.value, group = E.coli.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.coli.reshape.c5, aes(x=E.coli.reshape.c5$time, y = E.coli.reshape.c5$log10.exp.value, group = E.coli.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = E.coli.reshape.c6, aes(x=E.coli.reshape.c6$time, y = E.coli.reshape.c6$log10.exp.value, group = E.coli.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = E.coli.reshape.c7, aes(x=E.coli.reshape.c7$time, y = E.coli.reshape.c7$log10.exp.value, group = E.coli.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=5, newpage=T)

#S.mar.type - 12 clusters
q1 = ggplot(data = S.mar.type.reshape.c1, aes(x=S.mar.type.reshape.c1$time, y = S.mar.type.reshape.c1$log10.exp.value, group = S.mar.type.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = S.mar.type.reshape.c2, aes(x=S.mar.type.reshape.c2$time, y = S.mar.type.reshape.c2$log10.exp.value, group = S.mar.type.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = S.mar.type.reshape.c3, aes(x=S.mar.type.reshape.c3$time, y = S.mar.type.reshape.c3$log10.exp.value, group = S.mar.type.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = S.mar.type.reshape.c4, aes(x=S.mar.type.reshape.c4$time, y = S.mar.type.reshape.c4$log10.exp.value, group = S.mar.type.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = S.mar.type.reshape.c5, aes(x=S.mar.type.reshape.c5$time, y = S.mar.type.reshape.c5$log10.exp.value, group = S.mar.type.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = S.mar.type.reshape.c6, aes(x=S.mar.type.reshape.c6$time, y = S.mar.type.reshape.c6$log10.exp.value, group = S.mar.type.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = S.mar.type.reshape.c7, aes(x=S.mar.type.reshape.c7$time, y = S.mar.type.reshape.c7$log10.exp.value, group = S.mar.type.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = S.mar.type.reshape.c8, aes(x=S.mar.type.reshape.c8$time, y = S.mar.type.reshape.c8$log10.exp.value, group = S.mar.type.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = S.mar.type.reshape.c9, aes(x=S.mar.type.reshape.c9$time, y = S.mar.type.reshape.c9$log10.exp.value, group = S.mar.type.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = S.mar.type.reshape.c10, aes(x=S.mar.type.reshape.c10$time, y = S.mar.type.reshape.c10$log10.exp.value, group = S.mar.type.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = S.mar.type.reshape.c11, aes(x=S.mar.type.reshape.c11$time, y = S.mar.type.reshape.c11$log10.exp.value, group = S.mar.type.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = S.mar.type.reshape.c12, aes(x=S.mar.type.reshape.c12$time, y = S.mar.type.reshape.c12$log10.exp.value, group = S.mar.type.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#E.fae.live - 12 clusters
q1 = ggplot(data = E.fae.live.reshape.c1, aes(x=E.fae.live.reshape.c1$time, y = E.fae.live.reshape.c1$log10.exp.value, group = E.fae.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.fae.live.reshape.c2, aes(x=E.fae.live.reshape.c2$time, y = E.fae.live.reshape.c2$log10.exp.value, group = E.fae.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.fae.live.reshape.c3, aes(x=E.fae.live.reshape.c3$time, y = E.fae.live.reshape.c3$log10.exp.value, group = E.fae.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.fae.live.reshape.c4, aes(x=E.fae.live.reshape.c4$time, y = E.fae.live.reshape.c4$log10.exp.value, group = E.fae.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.fae.live.reshape.c5, aes(x=E.fae.live.reshape.c5$time, y = E.fae.live.reshape.c5$log10.exp.value, group = E.fae.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = E.fae.live.reshape.c6, aes(x=E.fae.live.reshape.c6$time, y = E.fae.live.reshape.c6$log10.exp.value, group = E.fae.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = E.fae.live.reshape.c7, aes(x=E.fae.live.reshape.c7$time, y = E.fae.live.reshape.c7$log10.exp.value, group = E.fae.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = E.fae.live.reshape.c8, aes(x=E.fae.live.reshape.c8$time, y = E.fae.live.reshape.c8$log10.exp.value, group = E.fae.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = E.fae.live.reshape.c9, aes(x=E.fae.live.reshape.c9$time, y = E.fae.live.reshape.c9$log10.exp.value, group = E.fae.live.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = E.fae.live.reshape.c10, aes(x=E.fae.live.reshape.c10$time, y = E.fae.live.reshape.c10$log10.exp.value, group = E.fae.live.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = E.fae.live.reshape.c11, aes(x=E.fae.live.reshape.c11$time, y = E.fae.live.reshape.c11$log10.exp.value, group = E.fae.live.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = E.fae.live.reshape.c12, aes(x=E.fae.live.reshape.c12$time, y = E.fae.live.reshape.c12$log10.exp.value, group = E.fae.live.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#P.rett.live - 15 clusters
q1 = ggplot(data = P.rett.live.reshape.c1, aes(x=P.rett.live.reshape.c1$time, y = P.rett.live.reshape.c1$log10.exp.value, group = P.rett.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = P.rett.live.reshape.c2, aes(x=P.rett.live.reshape.c2$time, y = P.rett.live.reshape.c2$log10.exp.value, group = P.rett.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = P.rett.live.reshape.c3, aes(x=P.rett.live.reshape.c3$time, y = P.rett.live.reshape.c3$log10.exp.value, group = P.rett.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = P.rett.live.reshape.c4, aes(x=P.rett.live.reshape.c4$time, y = P.rett.live.reshape.c4$log10.exp.value, group = P.rett.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = P.rett.live.reshape.c5, aes(x=P.rett.live.reshape.c5$time, y = P.rett.live.reshape.c5$log10.exp.value, group = P.rett.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = P.rett.live.reshape.c6, aes(x=P.rett.live.reshape.c6$time, y = P.rett.live.reshape.c6$log10.exp.value, group = P.rett.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = P.rett.live.reshape.c7, aes(x=P.rett.live.reshape.c7$time, y = P.rett.live.reshape.c7$log10.exp.value, group = P.rett.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = P.rett.live.reshape.c8, aes(x=P.rett.live.reshape.c8$time, y = P.rett.live.reshape.c8$log10.exp.value, group = P.rett.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = P.rett.live.reshape.c9, aes(x=P.rett.live.reshape.c9$time, y = P.rett.live.reshape.c9$log10.exp.value, group = P.rett.live.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = P.rett.live.reshape.c10, aes(x=P.rett.live.reshape.c10$time, y = P.rett.live.reshape.c10$log10.exp.value, group = P.rett.live.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = P.rett.live.reshape.c11, aes(x=P.rett.live.reshape.c11$time, y = P.rett.live.reshape.c11$log10.exp.value, group = P.rett.live.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = P.rett.live.reshape.c12, aes(x=P.rett.live.reshape.c12$time, y = P.rett.live.reshape.c12$log10.exp.value, group = P.rett.live.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q13 = ggplot(data = P.rett.live.reshape.c13, aes(x=P.rett.live.reshape.c13$time, y = P.rett.live.reshape.c13$log10.exp.value, group = P.rett.live.reshape.c13$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q14 = ggplot(data = P.rett.live.reshape.c14, aes(x=P.rett.live.reshape.c14$time, y = P.rett.live.reshape.c14$log10.exp.value, group = P.rett.live.reshape.c14$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q15 = ggplot(data = P.rett.live.reshape.c15, aes(x=P.rett.live.reshape.c15$time, y = P.rett.live.reshape.c15$log10.exp.value, group = P.rett.live.reshape.c15$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15, ncol=5, newpage=T)

#Ecc15 - 5 clusters
q1 = ggplot(data = Ecc15.reshape.c1, aes(x=Ecc15.reshape.c1$time, y = Ecc15.reshape.c1$log10.exp.value, group = Ecc15.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = Ecc15.reshape.c2, aes(x=Ecc15.reshape.c2$time, y = Ecc15.reshape.c2$log10.exp.value, group = Ecc15.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = Ecc15.reshape.c3, aes(x=Ecc15.reshape.c3$time, y = Ecc15.reshape.c3$log10.exp.value, group = Ecc15.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = Ecc15.reshape.c4, aes(x=Ecc15.reshape.c4$time, y = Ecc15.reshape.c4$log10.exp.value, group = Ecc15.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = Ecc15.reshape.c5, aes(x=Ecc15.reshape.c5$time, y = Ecc15.reshape.c5$log10.exp.value, group = Ecc15.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5, ncol=5, newpage=T)

#=========================================================================================

#(iv-2) Perform EBSeq-HMM on the expression of DE genes, excluding 0hr timepoint
M.luteus = M.luteus[,c(4:12)]; E.coli = E.coli[,(4:11)]; S.mar.type = S.mar.type[,c(4:12)]
E.fae.live = E.fae.live[,c(4:12)]; P.rett.live = P.rett.live[,c(4:12)]; Ecc15 = Ecc15[,c(4:12)]

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
for (j in 1:length(list.of.bacteria.name)){
    cat("EBSeq-HMM - We are working on the following comparison: 12hr-36hr-5.5d of ", as.character(list.of.bacteria.name[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 2){ #getting rid of Rep2 of E.coli 36hr
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
    
    write.table(de.list, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM/EBSeqHMM_DE_results_minus_0hr.tp_for_",list.of.bacteria.name[j],".txt", sep=""), quote=F)
}

#(v-2) Draw representation plots
count.table = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
library(ggplot2); library(grid); library(gridExtra)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM/")
for (s in 1:length(list.of.bacteria.name)){ 
    #For a given gene in the most_likely_path, identify the expression values in count.table and cbind
    most.likely.paths = read.table(paste("EBSeqHMM_DE_results_minus_0hr.tp_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
    #most.likely.paths = most.likely.paths[which(most.likely.paths$Max_PP > 0.55),] #getting rid of genes that have not so great categorization. Acutlly never used (5/23/2016)
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[s]))
    exp.list = c()
    for (n in 1:dim(most.likely.paths)[1]){
        gene.to.search = row.names(most.likely.paths[n,])
        gene.found = condition.count.table[which(row.names(condition.count.table) == gene.to.search),]
        exp.list = rbind(exp.list, gene.found)
    }
    if (s == 2){
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:5],1,mean),apply(exp.list[,6:8],1,mean))
    }
    else{
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:9],1,mean))
    }
    mean.bacteria.log10 = log10(mean.bacteria + 1) #This could be changed to centered mean
    most.likely.paths.bacteria = cbind(most.likely.paths,mean.bacteria.log10)
    colnames(most.likely.paths.bacteria) = c("paths","max_pp","12hr","36hr","5.5d")
    most.likely.paths.bacteria.reshape = reshape(most.likely.paths.bacteria, direction="long", varying = 3:5, idvar = "ids", ids = rownames(most.likely.paths.bacteria), timevar="time", v.names="log10.exp.value", times = c("12hr", "36hr", "5.5d"))
    most.likely.paths.bacteria = most.likely.paths.bacteria[order(most.likely.paths.bacteria$paths),]
    write.table(most.likely.paths.bacteria, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[s],"..minus.0hrtp.with.most.liekly.paths.txt",sep=""), quote=F, row.names =T, col.names=T)
    
    #Number of patterns
    pattern = levels(most.likely.paths.bacteria.reshape$paths)
    for (i in 1:length(pattern)){
        assign(paste0(list.of.bacteria.name[s],".reshape.c",i), most.likely.paths.bacteria.reshape[which(most.likely.paths.bacteria.reshape$paths == pattern[i]),])
    }
}

#Plot them here
#M.luteus - 8 clusters
q1 = ggplot(data = M.luteus.reshape.c1, aes(x=M.luteus.reshape.c1$time, y = M.luteus.reshape.c1$log10.exp.value, group = M.luteus.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = M.luteus.reshape.c2, aes(x=M.luteus.reshape.c2$time, y = M.luteus.reshape.c2$log10.exp.value, group = M.luteus.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = M.luteus.reshape.c3, aes(x=M.luteus.reshape.c3$time, y = M.luteus.reshape.c3$log10.exp.value, group = M.luteus.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = M.luteus.reshape.c4, aes(x=M.luteus.reshape.c4$time, y = M.luteus.reshape.c4$log10.exp.value, group = M.luteus.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = M.luteus.reshape.c5, aes(x=M.luteus.reshape.c5$time, y = M.luteus.reshape.c5$log10.exp.value, group = M.luteus.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = M.luteus.reshape.c6, aes(x=M.luteus.reshape.c6$time, y = M.luteus.reshape.c6$log10.exp.value, group = M.luteus.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = M.luteus.reshape.c7, aes(x=M.luteus.reshape.c7$time, y = M.luteus.reshape.c7$log10.exp.value, group = M.luteus.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = M.luteus.reshape.c8, aes(x=M.luteus.reshape.c8$time, y = M.luteus.reshape.c8$log10.exp.value, group = M.luteus.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=5, newpage=T)

#E.coli - 5 clusters
q1 = ggplot(data = E.coli.reshape.c1, aes(x=E.coli.reshape.c1$time, y = E.coli.reshape.c1$log10.exp.value, group = E.coli.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.coli.reshape.c2, aes(x=E.coli.reshape.c2$time, y = E.coli.reshape.c2$log10.exp.value, group = E.coli.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.coli.reshape.c3, aes(x=E.coli.reshape.c3$time, y = E.coli.reshape.c3$log10.exp.value, group = E.coli.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.coli.reshape.c4, aes(x=E.coli.reshape.c4$time, y = E.coli.reshape.c4$log10.exp.value, group = E.coli.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.coli.reshape.c5, aes(x=E.coli.reshape.c5$time, y = E.coli.reshape.c5$log10.exp.value, group = E.coli.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5, ncol=5, newpage=T)

#S.mar.type - 5 clusters
q1 = ggplot(data = S.mar.type.reshape.c1, aes(x=S.mar.type.reshape.c1$time, y = S.mar.type.reshape.c1$log10.exp.value, group = S.mar.type.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = S.mar.type.reshape.c2, aes(x=S.mar.type.reshape.c2$time, y = S.mar.type.reshape.c2$log10.exp.value, group = S.mar.type.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = S.mar.type.reshape.c3, aes(x=S.mar.type.reshape.c3$time, y = S.mar.type.reshape.c3$log10.exp.value, group = S.mar.type.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = S.mar.type.reshape.c4, aes(x=S.mar.type.reshape.c4$time, y = S.mar.type.reshape.c4$log10.exp.value, group = S.mar.type.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = S.mar.type.reshape.c5, aes(x=S.mar.type.reshape.c5$time, y = S.mar.type.reshape.c5$log10.exp.value, group = S.mar.type.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5, ncol=5, newpage=T)

#E.fae.live - 8 clusters
q1 = ggplot(data = E.fae.live.reshape.c1, aes(x=E.fae.live.reshape.c1$time, y = E.fae.live.reshape.c1$log10.exp.value, group = E.fae.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.fae.live.reshape.c2, aes(x=E.fae.live.reshape.c2$time, y = E.fae.live.reshape.c2$log10.exp.value, group = E.fae.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.fae.live.reshape.c3, aes(x=E.fae.live.reshape.c3$time, y = E.fae.live.reshape.c3$log10.exp.value, group = E.fae.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.fae.live.reshape.c4, aes(x=E.fae.live.reshape.c4$time, y = E.fae.live.reshape.c4$log10.exp.value, group = E.fae.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.fae.live.reshape.c5, aes(x=E.fae.live.reshape.c5$time, y = E.fae.live.reshape.c5$log10.exp.value, group = E.fae.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = E.fae.live.reshape.c6, aes(x=E.fae.live.reshape.c6$time, y = E.fae.live.reshape.c6$log10.exp.value, group = E.fae.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = E.fae.live.reshape.c7, aes(x=E.fae.live.reshape.c7$time, y = E.fae.live.reshape.c7$log10.exp.value, group = E.fae.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = E.fae.live.reshape.c8, aes(x=E.fae.live.reshape.c8$time, y = E.fae.live.reshape.c8$log10.exp.value, group = E.fae.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=5, newpage=T)

#P.rett.live - 8 clusters
q1 = ggplot(data = P.rett.live.reshape.c1, aes(x=P.rett.live.reshape.c1$time, y = P.rett.live.reshape.c1$log10.exp.value, group = P.rett.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = P.rett.live.reshape.c2, aes(x=P.rett.live.reshape.c2$time, y = P.rett.live.reshape.c2$log10.exp.value, group = P.rett.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = P.rett.live.reshape.c3, aes(x=P.rett.live.reshape.c3$time, y = P.rett.live.reshape.c3$log10.exp.value, group = P.rett.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = P.rett.live.reshape.c4, aes(x=P.rett.live.reshape.c4$time, y = P.rett.live.reshape.c4$log10.exp.value, group = P.rett.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = P.rett.live.reshape.c5, aes(x=P.rett.live.reshape.c5$time, y = P.rett.live.reshape.c5$log10.exp.value, group = P.rett.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = P.rett.live.reshape.c6, aes(x=P.rett.live.reshape.c6$time, y = P.rett.live.reshape.c6$log10.exp.value, group = P.rett.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = P.rett.live.reshape.c7, aes(x=P.rett.live.reshape.c7$time, y = P.rett.live.reshape.c7$log10.exp.value, group = P.rett.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = P.rett.live.reshape.c8, aes(x=P.rett.live.reshape.c8$time, y = P.rett.live.reshape.c8$log10.exp.value, group = P.rett.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=5, newpage=T)

#Ecc15 - 4 clusters
q1 = ggplot(data = Ecc15.reshape.c1, aes(x=Ecc15.reshape.c1$time, y = Ecc15.reshape.c1$log10.exp.value, group = Ecc15.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = Ecc15.reshape.c2, aes(x=Ecc15.reshape.c2$time, y = Ecc15.reshape.c2$log10.exp.value, group = Ecc15.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = Ecc15.reshape.c3, aes(x=Ecc15.reshape.c3$time, y = Ecc15.reshape.c3$log10.exp.value, group = Ecc15.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = Ecc15.reshape.c4, aes(x=Ecc15.reshape.c4$time, y = Ecc15.reshape.c4$log10.exp.value, group = Ecc15.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4, ncol=5, newpage=T)

#=========================================================================================
#(iii-3) Use variance stablized counts as the base to draw representation plots
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_all_time_points.txt", header=T); M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_all_time_points.txt", header=T); E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_all_time_points.txt", header=T); S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_all_time_points.txt", header=T); E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_all_time_points.txt", header=T); P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_all_time_points.txt", header=T); Ecc15.gene.only = Ecc15.gene.only[,1]

count.table.label = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/mega_RNA-seq_count_of_all_samples_Apr_2016_var_sta_transformed_May_2016.txt", header=T)
row.names(CountTable) = row.names(count.table.label)
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

#(iv-3) Draw representation plots using the existing EBSeq-HMM results
library(ggplot2); library(grid); library(gridExtra)
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM/")
for (s in 1:length(list.of.bacteria.name)){ 
    #For a given gene in the most_likely_path, identify the expression values in count.table and cbind
    most.likely.paths = read.table(paste("EBSeqHMM_DE_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
    #most.likely.paths = most.likely.paths[which(most.likely.paths$Max_PP > 0.55),] #getting rid of genes that have not so great categorization
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[s]))
    exp.list = c()
    for (n in 1:dim(most.likely.paths)[1]){
        gene.to.search = row.names(most.likely.paths[n,])
        gene.found = condition.count.table[which(row.names(condition.count.table) == gene.to.search),]
        exp.list = rbind(exp.list, gene.found)
    }
    if (s == 2){
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:8],1,mean),apply(exp.list[,9:11],1,mean))
    }
    else{
        mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:9],1,mean),apply(exp.list[,10:12],1,mean))
    }
    mean.bacteria.log10 = mean.bacteria
    most.likely.paths.bacteria = cbind(most.likely.paths,mean.bacteria.log10)
    colnames(most.likely.paths.bacteria) = c("paths","max_pp","0hr","12hr","36hr","5.5d")
    most.likely.paths.bacteria.reshape = reshape(most.likely.paths.bacteria, direction="long", varying = 3:6, idvar = "ids", ids = rownames(most.likely.paths.bacteria), timevar="time", v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
    most.likely.paths.bacteria = most.likely.paths.bacteria[order(most.likely.paths.bacteria$paths),]
    write.table(most.likely.paths.bacteria, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[s],".with.most.liekly.paths.based.on.VTS.txt",sep=""), quote=F, row.names =T, col.names=T)
    
    #Number of patterns
    pattern = levels(most.likely.paths.bacteria.reshape$paths)
    for (i in 1:length(pattern)){
        assign(paste0(list.of.bacteria.name[s],".reshape.c",i), most.likely.paths.bacteria.reshape[which(most.likely.paths.bacteria.reshape$paths == pattern[i]),])
    }
}

#Plot them here
#M.luteus - 12 clusters
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
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#E.coli - 7 clusters
q1 = ggplot(data = E.coli.reshape.c1, aes(x=E.coli.reshape.c1$time, y = E.coli.reshape.c1$log10.exp.value, group = E.coli.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.coli.reshape.c2, aes(x=E.coli.reshape.c2$time, y = E.coli.reshape.c2$log10.exp.value, group = E.coli.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.coli.reshape.c3, aes(x=E.coli.reshape.c3$time, y = E.coli.reshape.c3$log10.exp.value, group = E.coli.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.coli.reshape.c4, aes(x=E.coli.reshape.c4$time, y = E.coli.reshape.c4$log10.exp.value, group = E.coli.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.coli.reshape.c5, aes(x=E.coli.reshape.c5$time, y = E.coli.reshape.c5$log10.exp.value, group = E.coli.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = E.coli.reshape.c6, aes(x=E.coli.reshape.c6$time, y = E.coli.reshape.c6$log10.exp.value, group = E.coli.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = E.coli.reshape.c7, aes(x=E.coli.reshape.c7$time, y = E.coli.reshape.c7$log10.exp.value, group = E.coli.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=5, newpage=T)

#S.mar.type - 12 clusters
q1 = ggplot(data = S.mar.type.reshape.c1, aes(x=S.mar.type.reshape.c1$time, y = S.mar.type.reshape.c1$log10.exp.value, group = S.mar.type.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = S.mar.type.reshape.c2, aes(x=S.mar.type.reshape.c2$time, y = S.mar.type.reshape.c2$log10.exp.value, group = S.mar.type.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = S.mar.type.reshape.c3, aes(x=S.mar.type.reshape.c3$time, y = S.mar.type.reshape.c3$log10.exp.value, group = S.mar.type.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = S.mar.type.reshape.c4, aes(x=S.mar.type.reshape.c4$time, y = S.mar.type.reshape.c4$log10.exp.value, group = S.mar.type.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = S.mar.type.reshape.c5, aes(x=S.mar.type.reshape.c5$time, y = S.mar.type.reshape.c5$log10.exp.value, group = S.mar.type.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = S.mar.type.reshape.c6, aes(x=S.mar.type.reshape.c6$time, y = S.mar.type.reshape.c6$log10.exp.value, group = S.mar.type.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = S.mar.type.reshape.c7, aes(x=S.mar.type.reshape.c7$time, y = S.mar.type.reshape.c7$log10.exp.value, group = S.mar.type.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = S.mar.type.reshape.c8, aes(x=S.mar.type.reshape.c8$time, y = S.mar.type.reshape.c8$log10.exp.value, group = S.mar.type.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = S.mar.type.reshape.c9, aes(x=S.mar.type.reshape.c9$time, y = S.mar.type.reshape.c9$log10.exp.value, group = S.mar.type.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = S.mar.type.reshape.c10, aes(x=S.mar.type.reshape.c10$time, y = S.mar.type.reshape.c10$log10.exp.value, group = S.mar.type.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = S.mar.type.reshape.c11, aes(x=S.mar.type.reshape.c11$time, y = S.mar.type.reshape.c11$log10.exp.value, group = S.mar.type.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = S.mar.type.reshape.c12, aes(x=S.mar.type.reshape.c12$time, y = S.mar.type.reshape.c12$log10.exp.value, group = S.mar.type.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#E.fae.live - 12 clusters
q1 = ggplot(data = E.fae.live.reshape.c1, aes(x=E.fae.live.reshape.c1$time, y = E.fae.live.reshape.c1$log10.exp.value, group = E.fae.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = E.fae.live.reshape.c2, aes(x=E.fae.live.reshape.c2$time, y = E.fae.live.reshape.c2$log10.exp.value, group = E.fae.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = E.fae.live.reshape.c3, aes(x=E.fae.live.reshape.c3$time, y = E.fae.live.reshape.c3$log10.exp.value, group = E.fae.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = E.fae.live.reshape.c4, aes(x=E.fae.live.reshape.c4$time, y = E.fae.live.reshape.c4$log10.exp.value, group = E.fae.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = E.fae.live.reshape.c5, aes(x=E.fae.live.reshape.c5$time, y = E.fae.live.reshape.c5$log10.exp.value, group = E.fae.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = E.fae.live.reshape.c6, aes(x=E.fae.live.reshape.c6$time, y = E.fae.live.reshape.c6$log10.exp.value, group = E.fae.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = E.fae.live.reshape.c7, aes(x=E.fae.live.reshape.c7$time, y = E.fae.live.reshape.c7$log10.exp.value, group = E.fae.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = E.fae.live.reshape.c8, aes(x=E.fae.live.reshape.c8$time, y = E.fae.live.reshape.c8$log10.exp.value, group = E.fae.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = E.fae.live.reshape.c9, aes(x=E.fae.live.reshape.c9$time, y = E.fae.live.reshape.c9$log10.exp.value, group = E.fae.live.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = E.fae.live.reshape.c10, aes(x=E.fae.live.reshape.c10$time, y = E.fae.live.reshape.c10$log10.exp.value, group = E.fae.live.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = E.fae.live.reshape.c11, aes(x=E.fae.live.reshape.c11$time, y = E.fae.live.reshape.c11$log10.exp.value, group = E.fae.live.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = E.fae.live.reshape.c12, aes(x=E.fae.live.reshape.c12$time, y = E.fae.live.reshape.c12$log10.exp.value, group = E.fae.live.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12, ncol=5, newpage=T)

#P.rett.live - 15 clusters
q1 = ggplot(data = P.rett.live.reshape.c1, aes(x=P.rett.live.reshape.c1$time, y = P.rett.live.reshape.c1$log10.exp.value, group = P.rett.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = P.rett.live.reshape.c2, aes(x=P.rett.live.reshape.c2$time, y = P.rett.live.reshape.c2$log10.exp.value, group = P.rett.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = P.rett.live.reshape.c3, aes(x=P.rett.live.reshape.c3$time, y = P.rett.live.reshape.c3$log10.exp.value, group = P.rett.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = P.rett.live.reshape.c4, aes(x=P.rett.live.reshape.c4$time, y = P.rett.live.reshape.c4$log10.exp.value, group = P.rett.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = P.rett.live.reshape.c5, aes(x=P.rett.live.reshape.c5$time, y = P.rett.live.reshape.c5$log10.exp.value, group = P.rett.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q6 = ggplot(data = P.rett.live.reshape.c6, aes(x=P.rett.live.reshape.c6$time, y = P.rett.live.reshape.c6$log10.exp.value, group = P.rett.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q7 = ggplot(data = P.rett.live.reshape.c7, aes(x=P.rett.live.reshape.c7$time, y = P.rett.live.reshape.c7$log10.exp.value, group = P.rett.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q8 = ggplot(data = P.rett.live.reshape.c8, aes(x=P.rett.live.reshape.c8$time, y = P.rett.live.reshape.c8$log10.exp.value, group = P.rett.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q9 = ggplot(data = P.rett.live.reshape.c9, aes(x=P.rett.live.reshape.c9$time, y = P.rett.live.reshape.c9$log10.exp.value, group = P.rett.live.reshape.c9$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q10 = ggplot(data = P.rett.live.reshape.c10, aes(x=P.rett.live.reshape.c10$time, y = P.rett.live.reshape.c10$log10.exp.value, group = P.rett.live.reshape.c10$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q11 = ggplot(data = P.rett.live.reshape.c11, aes(x=P.rett.live.reshape.c11$time, y = P.rett.live.reshape.c11$log10.exp.value, group = P.rett.live.reshape.c11$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q12 = ggplot(data = P.rett.live.reshape.c12, aes(x=P.rett.live.reshape.c12$time, y = P.rett.live.reshape.c12$log10.exp.value, group = P.rett.live.reshape.c12$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q13 = ggplot(data = P.rett.live.reshape.c13, aes(x=P.rett.live.reshape.c13$time, y = P.rett.live.reshape.c13$log10.exp.value, group = P.rett.live.reshape.c13$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q14 = ggplot(data = P.rett.live.reshape.c14, aes(x=P.rett.live.reshape.c14$time, y = P.rett.live.reshape.c14$log10.exp.value, group = P.rett.live.reshape.c14$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q15 = ggplot(data = P.rett.live.reshape.c15, aes(x=P.rett.live.reshape.c15$time, y = P.rett.live.reshape.c15$log10.exp.value, group = P.rett.live.reshape.c15$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15, ncol=5, newpage=T)

#Ecc15 - 5 clusters
q1 = ggplot(data = Ecc15.reshape.c1, aes(x=Ecc15.reshape.c1$time, y = Ecc15.reshape.c1$log10.exp.value, group = Ecc15.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
q2 = ggplot(data = Ecc15.reshape.c2, aes(x=Ecc15.reshape.c2$time, y = Ecc15.reshape.c2$log10.exp.value, group = Ecc15.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q3 = ggplot(data = Ecc15.reshape.c3, aes(x=Ecc15.reshape.c3$time, y = Ecc15.reshape.c3$log10.exp.value, group = Ecc15.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q4 = ggplot(data = Ecc15.reshape.c4, aes(x=Ecc15.reshape.c4$time, y = Ecc15.reshape.c4$log10.exp.value, group = Ecc15.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
q5 = ggplot(data = Ecc15.reshape.c5, aes(x=Ecc15.reshape.c5$time, y = Ecc15.reshape.c5$log10.exp.value, group = Ecc15.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
grid.arrange(q1,q2,q3,q4,q5, ncol=5, newpage=T)


#=============
# par(mfrow=c(1,3))
# most.likely.paths.M.luteus.c1 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "Up-Up-Up"),]
# most.likely.paths.M.luteus.c1.exp.only = most.likely.paths.M.luteus.c1[,3:6]
# colnames(most.likely.paths.M.luteus.c1.exp.only) = c("0hr","12hr","36hr","5.5d")
# most.likely.paths.M.luteus.c1.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c1.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c1.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
# ggplot(data = most.likely.paths.M.luteus.c1.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c1.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c1.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c1.exp.only.reshaped$id)) + geom_line() + theme_classic()
# 
# most.likely.paths.M.luteus.c2 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "Up-Up-Down"),]
# most.likely.paths.M.luteus.c2.exp.only = most.likely.paths.M.luteus.c2[,3:6]
# colnames(most.likely.paths.M.luteus.c2.exp.only) = c("0hr","12hr","36hr","5.5d")
# most.likely.paths.M.luteus.c2.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c2.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c2.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
# ggplot(data = most.likely.paths.M.luteus.c2.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c2.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c2.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c2.exp.only.reshaped$id)) + geom_line() + theme_classic()
# 
# most.likely.paths.M.luteus.c3 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "-Down-Down"),]
# most.likely.paths.M.luteus.c3.exp.only = most.likely.paths.M.luteus.c3[,3:6]
# colnames(most.likely.paths.M.luteus.c3.exp.only) = c("0hr","12hr","36hr","5.5d")
# most.likely.paths.M.luteus.c3.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c3.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c3.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
# ggplot(data = most.likely.paths.M.luteus.c3.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c3.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c3.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c3.exp.only.reshaped$id)) + geom_line() + theme_classic()
# 


#(vi) Use the good old cluster package and plot the representation graph -- it's shit (5/23/2016)
# library(cluster) 
# for (t in 1:length(list.of.bacteria.name)){
#     # I. Partitioning
#     condition.count.table = list.of.bacteria[[t]] #Determine the number of clusters
# #     wss <- (nrow(condition.count.table)-1)*sum(apply(condition.count.table,2,var)) 
# #     for (i in 2:15) {wss[i] <- sum(kmeans(condition.count.table,centers=i)$withinss) }
# #     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/", list.of.bacteria.name[t], "_kmeans_WSS_based_on_VST.pdf", sep =""))
# #     plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") #the number plateaus at about 10
# #     dev.off()
#     
#     # II. K-Means Cluster Analysis
#     number.of.clusters = floor(sqrt(dim(condition.count.table)[1])/2) #equation from Yasir
#     fit <- kmeans(condition.count.table, number.of.clusters) # 10 cluster solution
#     aggregate(condition.count.table, by=list(fit$cluster),FUN=mean) # get cluster means 
#     condition.count.table <- cbind(condition.count.table, fit$cluster) # append cluster assignment
#     colnames(condition.count.table)[ncol(condition.count.table)] <- "cluster_info"
#     # vary parameters for most readable graph
#     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/", list.of.bacteria.name[t], "_kmeans_nK",number.of.clusters,"_cluster_based_on_VST.pdf", sep=""))
#     clusplot(condition.count.table, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
#     dev.off()
#     
#     # III. Ward Hierarchical Clustering
#     d <- dist(condition.count.table, method = "euclidean") # distance matrix
#     fit <- hclust(d, method="ward")
#     groups <- cutree(fit, k=number.of.clusters) # cut tree into number.of.clusters
#     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/", list.of.bacteria.name[t], "_kmeans_nK",number.of.clusters, "_hierarchical_cluster_based_on_VST.pdf", sep=""))
#     plot(fit) # display dendogram
#     rect.hclust(fit, k=number.of.clusters, border="red")# draw dendogram with red borders around the number.of.clusters
#     dev.off()
#     
#     write.table(condition.count.table, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/", list.of.bacteria.name[t], "_kmeans_nK",number.of.clusters, "_clustering_info_based_on_VST.txt", sep=""), quote=F, row.names = T, col.names = T)
# }
# 
# list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
# list.of.bacteria.file.name = c("M.luteus_kmeans_nK4_", "E.coli_kmeans_nK3_", "S.mar.type_kmeans_nK4_", "E.fae.live_kmeans_nK4_", "P.rett.live_kmeans_nK5_", "Ecc15_kmeans_nK2_")
# library(ggplot2); library(grid); library(gridExtra)
# for (s in 1:length(list.of.bacteria.name)){
#     bacteria.w.clustering.value = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/", list.of.bacteria.file.name[s], "clustering_info_based_on_VST.txt", sep=""), header=T)
#     if (s == 2){ #for E.coli missing one of the 36hr replicates
#         mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:5],1,mean),apply(bacteria.w.clustering.value[,6:8],1,mean),apply(bacteria.w.clustering.value[,9:11],1,mean),bacteria.w.clustering.value[,12])
#     }
#     else {
#         mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:6],1,mean),apply(bacteria.w.clustering.value[,7:9],1,mean),apply(bacteria.w.clustering.value[,10:12],1,mean),bacteria.w.clustering.value[,13])
#     }
#     colnames(mean.bacteria) = c("0hr","12hr","36hr","5.5d","cluster_number")
#     zero = data.frame(mean.bacteria[,1]); twelve = data.frame(mean.bacteria[,2]); thirty.six = data.frame(mean.bacteria[,3]); five.half = data.frame(mean.bacteria[,4])
#     zero.tag = cbind(zero, rep("0hr",dim(zero)[1]), mean.bacteria[,5]); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]), mean.bacteria[,5])
#     ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1]), mean.bacteria[,5]); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]), mean.bacteria[,5])
#     colnames(zero.tag) = c("log10Exp", "Time","Cluster");colnames(twelve.tag) = c("log10Exp", "Time","Cluster"); colnames(ts.tag) = c("log10Exp", "Time","Cluster"); colnames(fs.tag) = c("log10Exp", "time","Cluster")
#     total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.bacteria),4)); total = cbind(total,gene.names)
#     total$Time <- factor(total$Time); total$Cluster = factor(total$Cluster); total$gene.names = factor(total$gene.names)
#     
#     #Number of patterns
#     pattern = levels(total$Cluster)
#     for (i in 1:length(pattern)){
#         assign(paste0(list.of.bacteria.name[s],".reshape.c",i), total[which(total$Cluster == pattern[i]),])
#     }
# }
# 
# #M.luteus - 4 clusters
# q1 = ggplot(data = M.luteus.reshape.c1, aes(x=M.luteus.reshape.c1$Time, y = M.luteus.reshape.c1$log10Exp, group = M.luteus.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = M.luteus.reshape.c2, aes(x=M.luteus.reshape.c2$Time, y = M.luteus.reshape.c2$log10Exp, group = M.luteus.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q3 = ggplot(data = M.luteus.reshape.c3, aes(x=M.luteus.reshape.c3$Time, y = M.luteus.reshape.c3$log10Exp, group = M.luteus.reshape.c3$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q4 = ggplot(data = M.luteus.reshape.c4, aes(x=M.luteus.reshape.c4$Time, y = M.luteus.reshape.c4$log10Exp, group = M.luteus.reshape.c4$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2,q3,q4, ncol=5, newpage=T)
# 
# #E.coli - 3 clusters
# q1 = ggplot(data = E.coli.reshape.c1, aes(x=E.coli.reshape.c1$Time, y = E.coli.reshape.c1$log10Exp, group = E.coli.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = E.coli.reshape.c2, aes(x=E.coli.reshape.c2$Time, y = E.coli.reshape.c2$log10Exp, group = E.coli.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q3 = ggplot(data = E.coli.reshape.c3, aes(x=E.coli.reshape.c3$Time, y = E.coli.reshape.c3$log10Exp, group = E.coli.reshape.c3$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2,q3, ncol=5, newpage=T)
# 
# #S.mar.type - 4 clusters
# q1 = ggplot(data = S.mar.type.reshape.c1, aes(x=S.mar.type.reshape.c1$Time, y = S.mar.type.reshape.c1$log10Exp, group = S.mar.type.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = S.mar.type.reshape.c2, aes(x=S.mar.type.reshape.c2$Time, y = S.mar.type.reshape.c2$log10Exp, group = S.mar.type.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q3 = ggplot(data = S.mar.type.reshape.c3, aes(x=S.mar.type.reshape.c3$Time, y = S.mar.type.reshape.c3$log10Exp, group = S.mar.type.reshape.c3$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q4 = ggplot(data = S.mar.type.reshape.c4, aes(x=S.mar.type.reshape.c4$Time, y = S.mar.type.reshape.c4$log10Exp, group = S.mar.type.reshape.c4$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2,q3,q4, ncol=5, newpage=T)
# 
# #E.fae.live - 4 clusters
# q1 = ggplot(data = E.fae.live.reshape.c1, aes(x=E.fae.live.reshape.c1$Time, y = E.fae.live.reshape.c1$log10Exp, group = E.fae.live.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = E.fae.live.reshape.c2, aes(x=E.fae.live.reshape.c2$Time, y = E.fae.live.reshape.c2$log10Exp, group = E.fae.live.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q3 = ggplot(data = E.fae.live.reshape.c3, aes(x=E.fae.live.reshape.c3$Time, y = E.fae.live.reshape.c3$log10Exp, group = E.fae.live.reshape.c3$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q4 = ggplot(data = E.fae.live.reshape.c4, aes(x=E.fae.live.reshape.c4$Time, y = E.fae.live.reshape.c4$log10Exp, group = E.fae.live.reshape.c4$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2,q3,q4, ncol=5, newpage=T)
# 
# #P.rett.live - 5 clusters
# q1 = ggplot(data = P.rett.live.reshape.c1, aes(x=P.rett.live.reshape.c1$Time, y = P.rett.live.reshape.c1$log10Exp, group = P.rett.live.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = P.rett.live.reshape.c2, aes(x=P.rett.live.reshape.c2$Time, y = P.rett.live.reshape.c2$log10Exp, group = P.rett.live.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q3 = ggplot(data = P.rett.live.reshape.c3, aes(x=P.rett.live.reshape.c3$Time, y = P.rett.live.reshape.c3$log10Exp, group = P.rett.live.reshape.c3$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q4 = ggplot(data = P.rett.live.reshape.c4, aes(x=P.rett.live.reshape.c4$Time, y = P.rett.live.reshape.c4$log10Exp, group = P.rett.live.reshape.c4$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# q5 = ggplot(data = P.rett.live.reshape.c5, aes(x=P.rett.live.reshape.c5$Time, y = P.rett.live.reshape.c5$log10Exp, group = P.rett.live.reshape.c5$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2,q3,q4,q5, ncol=5, newpage=T)
# 
# #Ecc15 - 2 clusters
# q1 = ggplot(data = Ecc15.reshape.c1, aes(x=Ecc15.reshape.c1$Time, y = Ecc15.reshape.c1$log10Exp, group = Ecc15.reshape.c1$gene.names)) + geom_line() + theme_classic() +xlab("Time") + ylab("log10(expression)")
# q2 = ggplot(data = Ecc15.reshape.c2, aes(x=Ecc15.reshape.c2$Time, y = Ecc15.reshape.c2$log10Exp, group = Ecc15.reshape.c2$gene.names)) + geom_line() + theme_classic()+xlab("Time") + ylab("log10(expression)")
# grid.arrange(q1,q2, ncol=5, newpage=T)



