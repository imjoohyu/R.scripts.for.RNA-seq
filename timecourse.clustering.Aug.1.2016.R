#################################################################################################################
#Clustering analysis of RNA-seq data for kinetics/time-course
#8/1/2016
#Joo Hyun Im (ji72)

#Done in my laptop
rm(list=ls(all=TRUE))

#(i) Separate the table into 10 different conditions and pull out only the genes that are DE in at least one comparison
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
de.list = de.list[,c(1,2,9:52)]; M.luteus = de.list[,c(1:2,3:8)]; E.coli = de.list[,c(1:2,9:14)]; S.mar.type = de.list[,c(1:2,15:20)]; E.fae.live = de.list[,c(1:2,21:26)]; P.rett.live =de.list[,c(1:2,27:32)]; Ecc15 = de.list[,c(1:2,33:38)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

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


#(ii) Here's the list of DE genes that are significant in at least one comparison of time in each infection condition
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_at_least_one_time_point.txt", header=T); M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_at_least_one_time_point.txt", header=T); E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_at_least_one_time_point.txt", header=T); S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_at_least_one_time_point.txt", header=T); E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_at_least_one_time_point.txt", header=T); P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_at_least_one_time_point.txt", header=T); Ecc15.gene.only = Ecc15.gene.only[,1]

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

#(iii) Transform the data to log2(count+1)
M.luteus.transformed = log2(M.luteus + 1); E.coli.transformed = log2(E.coli + 1)
S.mar.type.transformed = log2(S.mar.type + 1); E.fae.live.transformed = log2(E.fae.live + 1)
P.rett.live.transformed = log2(P.rett.live + 1); Ecc15.transformed = log2(Ecc15 + 1)
list.of.bacteria = list(as.matrix(M.luteus.transformed), as.matrix(E.coli.transformed),as.matrix(S.mar.type.transformed), as.matrix(E.fae.live.transformed),as.matrix(P.rett.live.transformed), as.matrix(Ecc15.transformed))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")


#(iv) Cluster the genes based on expression using the cluster package
library(cluster) 
for (t in 1:length(list.of.bacteria.name)){
    # I. Partitioning
    condition.count.table = list.of.bacteria[[t]] #Determine the number of clusters
    number.of.cluster = round(sqrt(dim(condition.count.table)[1])/2)
    fit <- kmeans(condition.count.table, as.numeric(number.of.cluster)) # number of cluster determined by sqrt(# DEGs)/2
    aggregate(condition.count.table, by=list(fit$cluster),FUN=mean) # get cluster means 
    condition.count.table <- cbind(condition.count.table, fit$cluster) # append cluster assignment
    colnames(condition.count.table)[ncol(condition.count.table)] <- "cluster_info"
    # vary parameters for most readable graph
    #pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_cluster_based_on_VST.pdf", sep=""))
    #clusplot(condition.count.table, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
    #dev.off()
    
    # II. Ward Hierarchical Clustering
    d <- dist(condition.count.table, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward")
    groups <- cutree(fit, k=number.of.cluster) # cut tree into the number of cluster determined by sqrt(# DEGs)/2
    #pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_hierarchical_cluster_based_on_VST.pdf", sep=""))
    plot(fit) # display dendogram
    rect.hclust(fit, k=number.of.cluster, border="red")# draw dendogram with red borders around the X clusters
    #dev.off()
    
    # III. heatmap
    condition.count.table.scaled = as.matrix(scale(condition.count.table)) #scale data to mean=0, sd=1 and convert to matrix
    heatmap(condition.count.table.scaled, hclustfun = function(d) hclust(d,method="ward"))
        #write.table(condition.count.table, file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_clustering_info_based_on_VST.txt", sep=""), quote=F, row.names = T, col.names = T)
    
}






