#Network construction based on http://khughitt.github.io/2016-iscb-dc-rsg-workshop-presentation
#https://github.com/iscb-dc-rsg/2016-summer-workshop/tree/master/3B-Hughitt-RNASeq-Coex-Network-Analysis/tutorial
#May 17th, 2017
#Joo Hyun Im (ji72)

#delete any previous input
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")

library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library(WGCNA)
library('gplots')
library(RColorBrewer)


#1. Data pre-processing: select only 12hr samples
#low count genes have been filtered out, normalized
count_data = read.table("edgeR_normalized_counts_from_all_genes_with_all_UCs_no_filter_cpmF.txt", header=T)
#DEGs only
sig_genes_data = read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
#only select DEGs from the count data
count_data_sig_only = count_data[which(sig_genes_data[,1] %in% rownames(count_data)),]

#"Most of the methods developed for co-expression network analysis and network inference were written for use with microarray data, including WGCNA!" and thus we need to make our data look like microarray data by log2 transformation. " This will transform our discrete, over-dispersed counts to a more Poisson-like continuous distribution."
log2CPM = log2(count_data_sig_only[,c(2:103)] +1)
log2CPM = cbind(rownames(count_data_sig_only), count_data_sig_only[,1], log2CPM)
colnames(log2CPM)[1] = "gene_id"; colnames(log2CPM)[2] = "gene_name"
#Only select 12hr samples
#log2CPM_12h = log2CPM[,c(1,2,6,9,12,15,18,21,24,27,28,29,30,31,34,40,43,46,49,52,55,58,61,62,63,64,65,68,74,77,80,83,86,89,92,95,96,97,98,99,102)] #2589x41
#log2CPM_12h_count = log2CPM_12h[,c(3:41)]
#Only select 12hr P.rett, S. aureus, P. sneebia, Db11, and PE
log2CPM_12h = log2CPM[,c(1,2,21,27,28,29,30,55,61,62,63,64,89,95,96,97,98)] #2589x17
log2CPM_12h_count = log2CPM_12h[,c(3:17)]

#Does it look okay? I am not sure...
x = melt(as.matrix(log2CPM_12h_count))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + geom_density()

#Check sample correlation
heatmap.2(cor(log2CPM_12h_count), trace='none', main='Sample correlations (log2-transformed)')



#design = read.table("../mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)

#2. Adjacency matrix construction
#"First, we need to select a similarity measure to use when comparing gene expression profiles."
#Similarity measures: Pearson correlation, Spearman correlation, Bi-weight Midcorrelation, Euclidean distance, Mutual information

#"Similarity measure which combines elements from Pearson correlation and Euclidean distance."
#from Hughitt
cordist <- function(dat) {
    cor_matrix  <- cor(t(dat))
    
    dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
    dist_matrix <- log1p(dist_matrix)
    dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
    
    sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}
sim_matrix <- cordist(log2CPM_12h_count) #2589 x 2589

#Visualize
heatmap_indices <- sample(nrow(sim_matrix), 500)
heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

#Construct adjacency matrix by converting similarity matrix to an adjacency matrix
adj_matrix = adjacency.fromSimilarity(sim_matrix, power=12, type='signed')
adj_matrix = matrix(adj_matrix, nrow=nrow(adj_matrix))
#rownames(adj_matrix) = rownames(log2CPM_12h_count) #by gene ID
#colnames(adj_matrix) = rownames(log2CPM_12h_count)
rownames(adj_matrix) = log2CPM_12h[,2] #by gene name
colnames(adj_matrix) = log2CPM_12h[,2]

#Visualize
heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)

#3. Network Module Detection
#Detect co-expression modules in the network by using hierarchical clustering followed by branch-cutting

gene_tree =  hclust(as.dist(1 - adj_matrix), method="average")
module_labels = cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,deepSplit=TRUE) #2589
module_colors = labels2colors(module_labels) #2589


#export the data in an igraph form
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
    library('igraph')
    
    # Determine filename to use
    if (is.null(filename)) {
        filename='network.graphml'
    }
    
    # TODO 2015/04/09
    # Add option to rescale correlations for each module before applying
    # threshold (this is simpler than the previous approach of trying to
    # determine a different threshold for each module)
    #
    # Still, modules with very low correlations should be given somewhat
    # less priority than those with very high correlations.
    
    #module_colors <- unique(nodeAttrDataFrame$color)
    #module_genes <- which(nodeAttrDataFrame$color == color)
    #module_adjmat <- adj_mat[module_genes,]
    #num_genes <- length(module_genes)
    
    # Adjust threshold if needed to limit remaining edges
    max_edges <- max_edge_ratio * nrow(adj_mat)
    
    edge_to_total_ratio <- max_edges / length(adj_mat)
    edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
    
    # Also choose a minimum threshold to make sure that at least some edges
    # are left
    min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
    
    threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
    
    # Remove edges with weights lower than the cutoff
    adj_mat[abs(adj_mat) < threshold] <- 0
    
    # Drop any genes with no edges (TODO: Make optional)
    orphaned <- (colSums(adj_mat) == 0)
    adj_mat <- adj_mat[!orphaned, !orphaned]
    
    # Also remove annotation entries
    if (!is.null(nodeAttr)) {
        nodeAttr <- nodeAttr[!orphaned]
    }
    if (!is.null(nodeAttrDataFrame)) {
        nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
    }
    
    # Keep track of non-positive edges and rescale to range 0,1
    is_zero     <- adj_mat == 0
    is_negative <- adj_mat < 0
    
    adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
    adj_mat[is_zero] <- 0
    adj_mat[is_negative] <- -adj_mat[is_negative]
    
    if (verbose) {
        message(sprintf("Outputting matrix with %d nodes and %d edges", 
                        nrow(adj_mat), sum(adj_mat > 0)))
    }
    
    # Create a new graph and add vertices
    # Weighted graph
    if (weighted) {
        g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
    } else {
        adj_mat[adj_mat != 0] <- 1
        g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
    }
    
    # Add single node annotation from vector
    if (!is.null(nodeAttr)) {
        g <- set.vertex.attribute(g, "attr", value=nodeAttr)
    }
    
    # Add node one or more node annotations from a data frame
    if (!is.null(nodeAttrDataFrame)) {
        for (colname in colnames(nodeAttrDataFrame)) {
            g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
        }
    }
    
    edge_correlation_negative <- c()
    
    # neg_correlations[edge_list]
    edge_list <- get.edgelist(g)
    
    for (i in 1:nrow(edge_list)) {
        from <- edge_list[i, 1]    
        to   <- edge_list[i, 2]    
    }
    
    # Save graph to a file
    write.graph(g, filename, format='graphml')
    
    # return igraph
    return(g)
}

export_network_to_graphml(adj_matrix, filename='../network.genename.virulent.only.graphml',threshold=0.3, nodeAttrDataFrame=NULL)

