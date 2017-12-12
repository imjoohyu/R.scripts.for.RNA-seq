#Analysis on the triad (PE, Db11, Prett)
#April 19th, 2017
#Joo Hyun Im (ji72)


rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq")

#Get the data (don't worry about the FDR cut-off)
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64
expression_data_FC_only = expression_data[,c(1,2,seq(3,64,2))] #ignoring the p-value #33

#Pick genes whose average expression is 2x (or whatever number) higher in the triad compared to the rest. 
expression_data_FC_only$average_triad = rowMeans(expression_data_FC_only[,c(18,26,27)])
expression_data_FC_only$average_the_rest = rowMeans(expression_data_FC_only[,-c(1,2,18,26,27,34)])
expression_data_FC_only$diff_in_average = abs(expression_data_FC_only$average_triad - expression_data_FC_only$average_the_rest) #degree

#more than 2-fold difference
expression_data_FC_only_2FC_diff = expression_data_FC_only[which(expression_data_FC_only$diff_in_average > 1),] #167
expression_data_FC_only_2FC_diff_sorted = expression_data_FC_only_2FC_diff[order(expression_data_FC_only_2FC_diff$diff_in_average, decreasing = T),] #167 genes
write.table(expression_data_FC_only_2FC_diff_sorted[,c(1:2,34,35,36)], file="response_to_infection_damage_and_stress/genes_whose_expression_2-fold_higher_in_triad.txt", quote=F, row.names = F, col.names = T)


#Overlap between JNK (GO:0007254) and genes highly expressed in triad
JNK = read.table("gene_sets/JNK_cascade_genes_based_on_AmiGO.txt", header=F, sep='\t')
JNK[,1] = sapply(JNK[,1], function(x) {sub("FB:","",x)} )
JNK = unique(JNK)
overlap1 = intersect(expression_data_FC_only_2FC_diff_sorted[,1], JNK[,1]); length(overlap1) #1

#Overlap between JAK-STAT cascade (GO:0007259) and genes highly expressed in triad
JAKSTAT = read.table("gene_sets/JAK-STAT_cascade_genes_based_on_AmiGO.txt", header=F, sep='\t')
JAKSTAT[,1] = sapply(JAKSTAT[,1], function(x) {sub("FB:","",x)} )
JAKSTAT = unique(JAKSTAT)
overlap2 = intersect(expression_data_FC_only_2FC_diff_sorted[,1], JAKSTAT[,1]); length(overlap2) #1

#Overlap between Toll (GO:0008063) and genes highly expressed in triad
Toll = read.table("gene_sets/Toll_signaling_pathway_genes_based_on_AmiGO.txt", header=F, sep='\t')
Toll[,1] = sapply(Toll[,1], function(x) {sub("FB:","",x)} )
Toll = unique(Toll)
overlap3 = intersect(expression_data_FC_only_2FC_diff_sorted[,1], Toll[,1]); length(overlap3) #1






