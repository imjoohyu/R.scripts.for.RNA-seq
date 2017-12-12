#Testing whether all pathogenic infection (regardless of Gram type) activates the Toll pathway
#Date: September 23, 2015
#ji72

#Prediction: All will induce Toll except E.coli, Ecc15, P.rett.heatkilled

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_mega_RNA-seq_Sept_2015/")
data = read.table("Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_FCgene_ID-name_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
list.of.Toll.genes = read.csv("list.of.Toll.csv", header=F)
list.of.Toll.genes.vector = as.vector(t(list.of.Toll.genes))

#Pull the Toll genes out of the data
data.Toll <- matrix(NA, nrow=length(list.of.Toll.genes.vector), ncol=dim(data)[2])
for (i in 1:length(list.of.Toll.genes.vector)){
    gene.to.compare = list.of.Toll.genes.vector[i]
    checking.against.full.list <- data[which(data$gene_id == gene.to.compare),]
    #print(checking.against.full.list)
    if (dim(checking.against.full.list)[1] > 0){
        for (j in 1:64){
            data.Toll[i,j] <- as.character(checking.against.full.list[1,j])
        }
    }
}
data.Toll = na.omit(data.Toll); colnames(data.Toll) = colnames(data); data.Toll = data.frame(data.Toll) #make the table quotation-free.

#Q1. Which Toll genes have changed expression regardless of conditions at 12hr? (all "Y"s)
data.Toll.12hr = data.Toll[, c(1,2,9,10,15,16,21,22,27,28,33,34,39,40,45,46,47,48,49,50,51,52,53,54,59,60)] #select for 12hr time point conditions
data.Toll.12hr.all.on.gene_id <- vector(); data.Toll.12hr.all.on.gene_name <- vector()
length.of.table <- as.numeric(dim(data.Toll.12hr)[1]); width.of.table <- as.numeric(dim(data.Toll.12hr)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...24
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 26
        indicator <- isTRUE(data.Toll.12hr[k,m] == "Y") #Significant
        if (indicator == F) break #non-significant results
    }
    if (m == width.of.table && indicator == T){
        #print(data.Toll.12hr[k,1])
        data.Toll.12hr.all.on.gene_id <- append(data.Toll.12hr.all.on.gene_id, as.character(data.Toll.12hr[k,1]))   #only saving the gene ids that gave non-significant p-val in all of the columns
        data.Toll.12hr.all.on.gene_name <- append(data.Toll.12hr.all.on.gene_name, as.character(data.Toll.12hr[k,2]))
        }
}
data.Toll.12hr.all.on = data.frame(cbind(data.Toll.12hr.all.on.gene_id, data.Toll.12hr.all.on.gene_name)); data.Toll.12hr.all.on
#results: 7 genes

#Q2. Sort the conditions by Gram-positive and Gram-negative
data.Toll.12hr.by.Gram.type = data.Toll[, c(1,2,9,10,27,28,53,54,45,46,1,2,15,16,21,22,39,40,33,34,59,60,47,48,49,50,51,52)] 
data.Toll.12hr.by.Gram.type.tr = data.frame(t(data.Toll.12hr.by.Gram.type))
write.table(data.Toll.12hr.by.Gram.type.tr, file="Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_Toll_genes_only.txt", sep = "\t", quote=F, col.names=F)

data.Toll.12hr.by.Gram.type.Y.N.only = data.Toll[, c(10,28,54,46,16,22,40,34,60,48,50,52)] #Only get the Y/Ns for clustering
rownames(data.Toll.12hr.by.Gram.type.Y.N.only) = data.Toll[,2]

#Clustering and Heatmap
require(cluster)
require(gplots)
require(stringr)
for (i in 1:dim(data.Toll.12hr.by.Gram.type.Y.N.only)[2]){
    data.Toll.12hr.by.Gram.type.Y.N.only[,i]<- str_replace_all(data.Toll.12hr.by.Gram.type.Y.N.only[,i], "Y", 1)
    data.Toll.12hr.by.Gram.type.Y.N.only[,i]<- str_replace_all(data.Toll.12hr.by.Gram.type.Y.N.only[,i], "N", 0)
}
data.Toll.12hr.by.Gram.type.Y.N.only
data.Toll.12hr.by.Gram.type.Y.N.only = data.matrix(data.Toll.12hr.by.Gram.type.Y.N.only)
# pdf(file="Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_Toll_genes_only_cluster.pdf", height=5, width=5) 
# heatmap(data.Toll.12hr.by.Gram.type.Y.N.only, distfun = dist, col=redblue(256), density.info="none", 
#           trace="none", dendrogram="row", scale="row", cexCol=0.7, labRow=NA, Colv=NULL, ylab=rownames(data.Toll.12hr.by.Gram.type.Y.N.only))
# #hclustfun = hclust
# dev.off()

#Another way of clustering -- Ward Hierarchical Clustering
data.Toll.12hr.by.Gram.type.Y.N.only.tr <- t(data.Toll.12hr.by.Gram.type.Y.N.only)
pdf(file="Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_Toll_genes_only_cluster.pdf", height=5, width=5)
d <- dist(data.Toll.12hr.by.Gram.type.Y.N.only.tr, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
dev.off()

# #Q3. Which Toll genes have changed expression when exposed to any Gram-positive bacteria at 12hr? (any "Y"s)
# m.seq = c(4,10,16,24); indicator <- NULL; data.Toll.12hr.on.Gram.positive.gene_id = vector(); data.Toll.12hr.on.Gram.positive.gene_name = vector(); counter =0
# for (k in 1:length.of.table){ #1, 2, 3, ...24
#     for (m in m.seq){ #M.luteus, E.fae.live, S.aureus, E.fae.hk
#         indicator <- isTRUE(data.Toll.12hr[k,m] == "Y") #Significant
#         if (indicator == T) { counter= counter+1 }
#         if (indicator == F) break #non-significant results
#     }
#     if (counter > 0){
#         data.Toll.12hr.on.Gram.positive.gene_id <- append(data.Toll.12hr.on.Gram.positive.gene_id, as.character(data.Toll.12hr[k,1]))   #only saving the gene ids that gave non-significant p-val in all of the columns
#         data.Toll.12hr.on.Gram.positive.gene_name <- append(data.Toll.12hr.on.Gram.positive.gene_name, as.character(data.Toll.12hr[k,2]))
#     }
#     counter=0
# }
# data.Toll.12hr.on.Gram.positive = data.frame(cbind(data.Toll.12hr.on.Gram.positive.gene_id, data.Toll.12hr.on.Gram.positive.gene_name)); data.Toll.12hr.on.Gram.positive 
# #Results: 14 genes

# #Q4. Which Toll genes have changed expression when exposed to any Gram-negative bacteria at 12hr? (any "Y"s)
# m.seq = c(16,22,34,40,48,50,52,60); indicator <- NULL; data.Toll.12hr.on.Gram.negative.gene_id = vector(); data.Toll.12hr.on.Gram.negative.gene_name = vector(); counter =0
# for (k in 1:length.of.table){ #1, 2, 3, ...24
#     for (m in m.seq){ #E.coli, S.mar.type, P.rett.live, Ecc15, P.sneebia, S.mar.db11, PE, P.rett.hk
#         indicator <- isTRUE(data.Toll.12hr[k,m] == "Y") #Significant
#         if (indicator == T) { counter= counter+1 }; print(counter)
#         if (indicator == F)  break #non-significant results
#     }
#     if (counter > 0){
#         data.Toll.12hr.on.Gram.negative.gene_id <- append(data.Toll.12hr.on.Gram.negative.gene_id, as.character(data.Toll.12hr[k,1]))   #only saving the gene ids that gave non-significant p-val in all of the columns
#         data.Toll.12hr.on.Gram.negative.gene_name <- append(data.Toll.12hr.on.Gram.negative.gene_name, as.character(data.Toll.12hr[k,2]))
#     }
#     counter=0
# }
# data.Toll.12hr.on.Gram.negative = data.frame(cbind(data.Toll.12hr.on.Gram.negative.gene_id, data.Toll.12hr.on.Gram.negative.gene_name)); data.Toll.12hr.on.Gram.negative
# #Results: 19 genes

#Q5. Which Toll genes have changed expression 
