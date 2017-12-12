#################################################################################################################
#Use EBSeq-HMM to make clusters of DE genes
#Updated on August 23, 2016
#Joo Hyun Im (ji72)
#################################################################################################################

#delete any previous input
rm(list=ls(all=TRUE))

#0. Separate the table into 10 different conditions and pull out only the genes that are DE in at least one comparison
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
de.list = de.list[,c(1,2,9:52)]; M.luteus = de.list[,c(1:2,3:8)]; E.coli = de.list[,c(1:2,9:14)]; S.mar.type = de.list[,c(1:2,15:20)]; E.fae.live = de.list[,c(1:2,21:26)]; P.rett.live =de.list[,c(1:2,27:32)]; Ecc15 = de.list[,c(1:2,33:38)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
#de.list below is from the original file, not the subset de.list[,c(1,2,9:52)]
E.fae.heatkilled = de.list[,c(1,2,53:58)]; P.rett.heatkilled = de.list[,c(1,2,59:64)]
list.of.bacteria = list(as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("E.fae.heatkilled", "P.rett.heatkilled")

for (i in 1:length(list.of.bacteria.name)){
    condition.count.table = list.of.bacteria[[i]]
    sig.DE.only = c()
    for (j in 1:dim(condition.count.table)[1]){
        number.of.positives = length(as.integer(grep("Y", condition.count.table[j,])))
        if (number.of.positives > 0){
            sig.DE.only = rbind(sig.DE.only, condition.count.table[j,])
        }
    }
    write.table(sig.DE.only, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[i],"_sig_DE_in_at_least_one_time_point.txt",sep=""), quote=F, row.names = F, col.names = T)
}


#1. Read in the data table of DE genes that are significant in at least one comparison of time point in each infection condition.
#For instance, if a gene was DE in 0hr-12hr OR 0hr-36hr OR 0hr-5.5d comparison, the gene was considered significant.
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_at_least_one_time_point.txt", header=T); M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_at_least_one_time_point.txt", header=T); E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_at_least_one_time_point.txt", header=T); S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_at_least_one_time_point.txt", header=T); E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_at_least_one_time_point.txt", header=T); P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_at_least_one_time_point.txt", header=T); Ecc15.gene.only = Ecc15.gene.only[,1]
E.fae.heatkilled.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.heatkilled_sig_DE_in_at_least_one_time_point.txt", header=T); E.fae.heatkilled.gene.only = E.fae.heatkilled.gene.only[,1]
P.rett.heatkilled.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.heatkilled_sig_DE_in_at_least_one_time_point.txt", header=T); P.rett.heatkilled.gene.only = P.rett.heatkilled.gene.only[,1]


#2. Pull out the expression values of DE genes from the counttable and match them with the gene names from #1.
#On August 23, 2016: I changed this to normalized count data because the EBSeq-HMM assumes that the library sizes are equal among samples.
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_with_all_UCs.txt", header=T) #Normalized count data (11911 x 103). Filtered.
CountTable = CountTable[,c(2:103)] #getting rid of the gene name
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
E.fae.heatkilled = cbind(unchallenged.average,E.fae.heatkilled.12hr,E.fae.heatkilled.36hr,E.fae.heatkilled.5.5d); E.fae.heatkilled = E.fae.heatkilled[row.names(E.fae.heatkilled) %in% E.fae.heatkilled.gene.only,]
P.rett.heatkilled = cbind(unchallenged.average,P.rett.heatkilled.12hr,P.rett.heatkilled.36hr,P.rett.heatkilled.5.5d); P.rett.heatkilled = P.rett.heatkilled[row.names(P.rett.heatkilled) %in% P.rett.heatkilled.gene.only,]


#3. Perform EBSeq-HMM on the raw expression data of DE genes
library(EBSeqHMM)
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

for (j in 1:length(list.of.bacteria.name)){
    cat("EBSeq-HMM - We are working on the following comparison: 0hr-12hr-36hr-5.5d of ", as.character(list.of.bacteria.name[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 2){ #getting rid of Rep2 of E.coli 36hr
        CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d")} #removing 36hr rep2 sample 
    else {
        CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d") 
    }    
    Conditions <- factor(CondVector, levels=c("t.0hr","t.12hr","t.36hr","t.5.5d"))
    
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=count.table, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
    EBSeqHMMGeneOut$MgAllMAPChar #classifies the pattern of gene expression over time per gene
    AllPaths = GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=TRUE) #Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
    gene.list.for.a.given.path = GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, cutoff=0.5, OnlyDynamic=TRUE,Paths=NULL) #Obtain confident gene calls for classifying genes into expression paths, PP>=0.5 is most likely path
    number.of.genes.in.each.category = gene.list.for.a.given.path$NumEach; number.of.genes.in.each.category #already sorted
    genes.in.each.category = gene.list.for.a.given.path$EachPathNames; genes.in.each.category
    de.list = GetDECalls(EBSeqHMMGeneOut,FDR=.05) #Another way to obtain DE gene/isoform list at a certain FDR (but not PPs)
    
    #For a specific gene, draw a plot of expression
    #PlotExp(count.table, Conditions, "Gene_1")
    
    write.table(de.list, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_",list.of.bacteria.name[j],".txt", sep=""), quote=F)
}


#4. Process EBSeq-HMM results and draw representation plots
#On 8/23, the data set was changed to use the normalized count data
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/")

#Determine the cutoff for the most likely path for each infection condition (added on 8/22/2016)
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

list.of.M.luteus.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_M.luteus.txt", header=T)
list.of.E.coli.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_E.coli.txt", header=T)
list.of.S.mar.type.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_S.mar.type.txt", header=T)
list.of.E.fae.live.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_E.fae.live.txt", header=T)
list.of.P.rett.live.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_P.rett.live.txt", header=T)
list.of.Ecc15.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_Ecc15.txt", header=T)
list.of.E.fae.heatkilled.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_E.fae.heatkilled.txt", header=T)
list.of.P.rett.heatkilled.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_normalized_count_results_for_P.rett.heatkilled.txt", header=T)
list.of.genes = list(as.matrix(list.of.M.luteus.genes), as.matrix(list.of.E.coli.genes), as.matrix(list.of.S.mar.type.genes), 
                     as.matrix(list.of.E.fae.live.genes), as.matrix(list.of.P.rett.live.genes), as.matrix(list.of.Ecc15.genes), as.matrix(list.of.E.fae.heatkilled.genes), as.matrix(list.of.P.rett.heatkilled.genes))

for (g in 1:length(list.of.genes)){
    list.of.genes.actual.results = c()
    sample = as.matrix(data.frame(list.of.bacteria[g]))
    sample.gene = as.matrix(data.frame(list.of.genes[g]))
    for (h in 1:dim(sample.gene)[1]){
        a = as.matrix(sample[which(row.names(sample) == row.names(sample.gene)[h]),])
        
        if (g == 2){ #E.coli
            one = median(a[c(1:3),1]); two = median(a[c(4:6),1]); three = median(a[c(7:8),1]); four = median(a[c(9:11),1])
        }
        else{ #Everything else
            one = median(a[c(1:3),1]); two = median(a[c(4:6),1]); three = median(a[c(7:9),1]); four = median(a[c(10:12),1])
        }
        
        if (one-two > 0){
            dir1="Down"
        }
        else {
            dir1="Up"
        }
        if (two-three >0){
            dir2="Down"
        }
        else{
            dir2="Up"
        }
        if (three-four > 0){
            dir3="Down"
        }
        else{
            dir3="Up"
        }
        full.direction = paste(dir1,"-",dir2,"-",dir3,sep="")
        list.of.genes.actual.results[h] = full.direction
    }
    final = cbind(sample.gene, list.of.genes.actual.results)
    yes.no=c()
    for (s in 1:dim(final)[1]){
        if(final[s,1] == final[s,3]){
            yes.no=rbind(yes.no,"Y")
        }
        else{
            yes.no=rbind(yes.no,"N")
        }
    }
    final =cbind (final, yes.no)
    write.table(final, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/EBSeqHMM_DE_results_checking_cutoff_for_", list.of.bacteria.name[g], "_normalized_data.txt", sep=""), quote=F)
}

# #Based on this comparison between EBSeq predictions and my simple math predictions, I determined the cutoff for each condition -- no longer used
# cut.off = as.numeric(c("0.5012", "0.5001", "0.5957", "0.5539", "0.5442","0.5373","0.6768","0.5463"))

#Filter out those below the cut.off probability and put together expression values with gene name and cluster types
number.of.cluster = c()
for (s in 1:length(list.of.bacteria.name)){ 
    #For a given gene in the most_likely_path, identify the expression values in count.table and cbind
    most.likely.paths = read.table(paste("EBSeqHMM_DE_normalized_count_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
    most.likely.paths = most.likely.paths[which(most.likely.paths$Max_PP > 0.5),] #getting rid of genes that have not so great categorization
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[s]))
    exp.list = c()
    for (n in 1:dim(most.likely.paths)[1]){
        gene.to.search = row.names(most.likely.paths[n,])
        gene.found = condition.count.table[which(row.names(condition.count.table) == gene.to.search),]
        exp.list = rbind(exp.list, gene.found)
    }
    if (s == 2){
        median.bacteria = cbind(apply(exp.list[,1:3],1,median),apply(exp.list[,4:6],1,median),apply(exp.list[,7:8],1,median),apply(exp.list[,9:11],1,median))
    }
    else{
        median.bacteria = cbind(apply(exp.list[,1:3],1,median),apply(exp.list[,4:6],1,median),apply(exp.list[,7:9],1,median),apply(exp.list[,10:12],1,median))
    }
    median.bacteria.log2 = log2(median.bacteria + 1) #transform data for plotting
    most.likely.paths.bacteria = cbind(most.likely.paths,median.bacteria.log2)
    colnames(most.likely.paths.bacteria) = c("paths","max_pp","0hr","12hr","36hr","5.5d")
    most.likely.paths.bacteria.reshape = reshape(most.likely.paths.bacteria, direction="long", varying = 3:6, idvar = "ids", ids = rownames(most.likely.paths.bacteria), timevar="time", v.names="log2.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
    most.likely.paths.bacteria = most.likely.paths.bacteria[order(most.likely.paths.bacteria$paths),]
    #write.table(most.likely.paths.bacteria, file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/",list.of.bacteria.name[s],".with.most.liekly.paths.based.on.normalized.counts.txt",sep=""), quote=F, row.names =T, col.names=T)
    
    #Number of patterns
    most.likely.paths.bacteria.reshape$paths = factor(most.likely.paths.bacteria.reshape$paths)
    pattern = levels(most.likely.paths.bacteria.reshape$paths)
    for (i in 1:length(pattern)){
        assign(paste0(list.of.bacteria.name[s],".reshape.c",i), most.likely.paths.bacteria.reshape[which(most.likely.paths.bacteria.reshape$paths == pattern[i]),])
    }
    
    cat("Number of clusters for ", list.of.bacteria.name[s], " is ", length(pattern), "\n")
    number.of.cluster[s] = as.integer(length(pattern))
}

#5. Plot the clusters
library(ggplot2); library(grid); library(gridExtra)

#M.luteus - 8 clusters
q1 = ggplot(data = M.luteus.reshape.c1, aes(x=M.luteus.reshape.c1$time, y = M.luteus.reshape.c1$log2.exp.value, group = M.luteus.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = M.luteus.reshape.c2, aes(x=M.luteus.reshape.c2$time, y = M.luteus.reshape.c2$log2.exp.value, group = M.luteus.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = M.luteus.reshape.c3, aes(x=M.luteus.reshape.c3$time, y = M.luteus.reshape.c3$log2.exp.value, group = M.luteus.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = M.luteus.reshape.c4, aes(x=M.luteus.reshape.c4$time, y = M.luteus.reshape.c4$log2.exp.value, group = M.luteus.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = M.luteus.reshape.c5, aes(x=M.luteus.reshape.c5$time, y = M.luteus.reshape.c5$log2.exp.value, group = M.luteus.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = M.luteus.reshape.c6, aes(x=M.luteus.reshape.c6$time, y = M.luteus.reshape.c6$log2.exp.value, group = M.luteus.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = M.luteus.reshape.c7, aes(x=M.luteus.reshape.c7$time, y = M.luteus.reshape.c7$log2.exp.value, group = M.luteus.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q8 = ggplot(data = M.luteus.reshape.c8, aes(x=M.luteus.reshape.c8$time, y = M.luteus.reshape.c8$log2.exp.value, group = M.luteus.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
print(grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=4, newpage=T))

#E.coli - 7 clusters
q1 = ggplot(data = E.coli.reshape.c1, aes(x=E.coli.reshape.c1$time, y = E.coli.reshape.c1$log2.exp.value, group = E.coli.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = E.coli.reshape.c2, aes(x=E.coli.reshape.c2$time, y = E.coli.reshape.c2$log2.exp.value, group = E.coli.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = E.coli.reshape.c3, aes(x=E.coli.reshape.c3$time, y = E.coli.reshape.c3$log2.exp.value, group = E.coli.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = E.coli.reshape.c4, aes(x=E.coli.reshape.c4$time, y = E.coli.reshape.c4$log2.exp.value, group = E.coli.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = E.coli.reshape.c5, aes(x=E.coli.reshape.c5$time, y = E.coli.reshape.c5$log2.exp.value, group = E.coli.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = E.coli.reshape.c6, aes(x=E.coli.reshape.c6$time, y = E.coli.reshape.c6$log2.exp.value, group = E.coli.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = E.coli.reshape.c7, aes(x=E.coli.reshape.c7$time, y = E.coli.reshape.c7$log2.exp.value, group = E.coli.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=4, newpage=T)

#S.mar.type - 7 clusters
q1 = ggplot(data = S.mar.type.reshape.c1, aes(x=S.mar.type.reshape.c1$time, y = S.mar.type.reshape.c1$log2.exp.value, group = S.mar.type.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = S.mar.type.reshape.c2, aes(x=S.mar.type.reshape.c2$time, y = S.mar.type.reshape.c2$log2.exp.value, group = S.mar.type.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = S.mar.type.reshape.c3, aes(x=S.mar.type.reshape.c3$time, y = S.mar.type.reshape.c3$log2.exp.value, group = S.mar.type.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = S.mar.type.reshape.c4, aes(x=S.mar.type.reshape.c4$time, y = S.mar.type.reshape.c4$log2.exp.value, group = S.mar.type.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = S.mar.type.reshape.c5, aes(x=S.mar.type.reshape.c5$time, y = S.mar.type.reshape.c5$log2.exp.value, group = S.mar.type.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = S.mar.type.reshape.c6, aes(x=S.mar.type.reshape.c6$time, y = S.mar.type.reshape.c6$log2.exp.value, group = S.mar.type.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = S.mar.type.reshape.c7, aes(x=S.mar.type.reshape.c7$time, y = S.mar.type.reshape.c7$log2.exp.value, group = S.mar.type.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=4, newpage=T)

#E.fae.live - 8 clusters
q1 = ggplot(data = E.fae.live.reshape.c1, aes(x=E.fae.live.reshape.c1$time, y = E.fae.live.reshape.c1$log2.exp.value, group = E.fae.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = E.fae.live.reshape.c2, aes(x=E.fae.live.reshape.c2$time, y = E.fae.live.reshape.c2$log2.exp.value, group = E.fae.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = E.fae.live.reshape.c3, aes(x=E.fae.live.reshape.c3$time, y = E.fae.live.reshape.c3$log2.exp.value, group = E.fae.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = E.fae.live.reshape.c4, aes(x=E.fae.live.reshape.c4$time, y = E.fae.live.reshape.c4$log2.exp.value, group = E.fae.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = E.fae.live.reshape.c5, aes(x=E.fae.live.reshape.c5$time, y = E.fae.live.reshape.c5$log2.exp.value, group = E.fae.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = E.fae.live.reshape.c6, aes(x=E.fae.live.reshape.c6$time, y = E.fae.live.reshape.c6$log2.exp.value, group = E.fae.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = E.fae.live.reshape.c7, aes(x=E.fae.live.reshape.c7$time, y = E.fae.live.reshape.c7$log2.exp.value, group = E.fae.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = E.fae.live.reshape.c8, aes(x=E.fae.live.reshape.c8$time, y = E.fae.live.reshape.c8$log2.exp.value, group = E.fae.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=4, newpage=T)

#P.rett.live - 8 clusters
q1 = ggplot(data = P.rett.live.reshape.c1, aes(x=P.rett.live.reshape.c1$time, y = P.rett.live.reshape.c1$log2.exp.value, group = P.rett.live.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = P.rett.live.reshape.c2, aes(x=P.rett.live.reshape.c2$time, y = P.rett.live.reshape.c2$log2.exp.value, group = P.rett.live.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = P.rett.live.reshape.c3, aes(x=P.rett.live.reshape.c3$time, y = P.rett.live.reshape.c3$log2.exp.value, group = P.rett.live.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = P.rett.live.reshape.c4, aes(x=P.rett.live.reshape.c4$time, y = P.rett.live.reshape.c4$log2.exp.value, group = P.rett.live.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = P.rett.live.reshape.c5, aes(x=P.rett.live.reshape.c5$time, y = P.rett.live.reshape.c5$log2.exp.value, group = P.rett.live.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = P.rett.live.reshape.c6, aes(x=P.rett.live.reshape.c6$time, y = P.rett.live.reshape.c6$log2.exp.value, group = P.rett.live.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = P.rett.live.reshape.c7, aes(x=P.rett.live.reshape.c7$time, y = P.rett.live.reshape.c7$log2.exp.value, group = P.rett.live.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q8 = ggplot(data = P.rett.live.reshape.c8, aes(x=P.rett.live.reshape.c8$time, y = P.rett.live.reshape.c8$log2.exp.value, group = P.rett.live.reshape.c8$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8, ncol=4, newpage=T)

#Ecc15 - 7 clusters
q1 = ggplot(data = Ecc15.reshape.c1, aes(x=Ecc15.reshape.c1$time, y = Ecc15.reshape.c1$log2.exp.value, group = Ecc15.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = Ecc15.reshape.c2, aes(x=Ecc15.reshape.c2$time, y = Ecc15.reshape.c2$log2.exp.value, group = Ecc15.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = Ecc15.reshape.c3, aes(x=Ecc15.reshape.c3$time, y = Ecc15.reshape.c3$log2.exp.value, group = Ecc15.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = Ecc15.reshape.c4, aes(x=Ecc15.reshape.c4$time, y = Ecc15.reshape.c4$log2.exp.value, group = Ecc15.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = Ecc15.reshape.c5, aes(x=Ecc15.reshape.c5$time, y = Ecc15.reshape.c5$log2.exp.value, group = Ecc15.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = Ecc15.reshape.c6, aes(x=Ecc15.reshape.c6$time, y = Ecc15.reshape.c6$log2.exp.value, group = Ecc15.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = Ecc15.reshape.c7, aes(x=Ecc15.reshape.c7$time, y = Ecc15.reshape.c7$log2.exp.value, group = Ecc15.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=4, newpage=T)

#E.faecalis heatkilled - 7 clusters
q1 = ggplot(data = E.fae.heatkilled.reshape.c1, aes(x=E.fae.heatkilled.reshape.c1$time, y = E.fae.heatkilled.reshape.c1$log2.exp.value, group = E.fae.heatkilled.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = E.fae.heatkilled.reshape.c2, aes(x=E.fae.heatkilled.reshape.c2$time, y = E.fae.heatkilled.reshape.c2$log2.exp.value, group = E.fae.heatkilled.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = E.fae.heatkilled.reshape.c3, aes(x=E.fae.heatkilled.reshape.c3$time, y = E.fae.heatkilled.reshape.c3$log2.exp.value, group = E.fae.heatkilled.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = E.fae.heatkilled.reshape.c4, aes(x=E.fae.heatkilled.reshape.c4$time, y = E.fae.heatkilled.reshape.c4$log2.exp.value, group = E.fae.heatkilled.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = E.fae.heatkilled.reshape.c5, aes(x=E.fae.heatkilled.reshape.c5$time, y = E.fae.heatkilled.reshape.c5$log2.exp.value, group = E.fae.heatkilled.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q6 = ggplot(data = E.fae.heatkilled.reshape.c6, aes(x=E.fae.heatkilled.reshape.c6$time, y = E.fae.heatkilled.reshape.c6$log2.exp.value, group = E.fae.heatkilled.reshape.c6$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q7 = ggplot(data = E.fae.heatkilled.reshape.c7, aes(x=E.fae.heatkilled.reshape.c7$time, y = E.fae.heatkilled.reshape.c7$log2.exp.value, group = E.fae.heatkilled.reshape.c7$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
grid.arrange(q1,q2,q3,q4,q5,q6,q7, ncol=4, newpage=T)

#P.rettgeri heatkilled - 5 clusters
q1 = ggplot(data = P.rett.heatkilled.reshape.c1, aes(x=P.rett.heatkilled.reshape.c1$time, y = P.rett.heatkilled.reshape.c1$log2.exp.value, group = P.rett.heatkilled.reshape.c1$ids)) + geom_line() + theme_classic() +xlab("Time") + ylab("log2(expression)")
q2 = ggplot(data = P.rett.heatkilled.reshape.c2, aes(x=P.rett.heatkilled.reshape.c2$time, y = P.rett.heatkilled.reshape.c2$log2.exp.value, group = P.rett.heatkilled.reshape.c2$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q3 = ggplot(data = P.rett.heatkilled.reshape.c3, aes(x=P.rett.heatkilled.reshape.c3$time, y = P.rett.heatkilled.reshape.c3$log2.exp.value, group = P.rett.heatkilled.reshape.c3$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q4 = ggplot(data = P.rett.heatkilled.reshape.c4, aes(x=P.rett.heatkilled.reshape.c4$time, y = P.rett.heatkilled.reshape.c4$log2.exp.value, group = P.rett.heatkilled.reshape.c4$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
q5 = ggplot(data = P.rett.heatkilled.reshape.c5, aes(x=P.rett.heatkilled.reshape.c5$time, y = P.rett.heatkilled.reshape.c5$log2.exp.value, group = P.rett.heatkilled.reshape.c5$ids)) + geom_line() + theme_classic()+xlab("Time") + ylab("log2(expression)")
print(grid.arrange(q1,q2,q3,q4,q5, ncol=4, newpage=T))
