#PE specificity
#March 31st, 2017

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq")
PE = read.table("specific.comparisons/genes_uniquely_regulated_by_each_condition/unique_genes_ignoring_time/uniquely_Upregulated_by_P.entomophila_ignoring_time.txt", header=T)

data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header =T) #2589 x 64

#Get the expression data of the PE-only genes
PE_data = data[which(data$gene_id %in% PE[,1] ),] #122
PE_data = droplevels(PE_data)
PE_data_sorted = PE_data[order(PE_data$log2FC.P.ento.12hr, decreasing=T),]

#And pick the ones that have very different DE patterns in PE than in other conditions
#Most of these are lowly expressed genes