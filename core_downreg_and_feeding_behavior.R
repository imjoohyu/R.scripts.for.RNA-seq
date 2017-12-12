#compare core downregulated genes with feeding behavior genes
#April 4, 2017
#Joo Hyun Im

core_downreg = read.table("finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
feeding = read.table("gene_sets/feeding_behavior_genes_based_on_AmiGO.txt", header=F, sep="\t")
feeding[,1] = sapply(feeding[,1], function(x) {sub("FB:","",x)} )
feeding = unique(feeding)

overlap = intersect(core_downreg[,1], feeding[,1])
length(overlap) #1

expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig_direction.txt", header=T) #FC, not count data, 11911 x 64
expression_data_subset = expression_data[match(feeding[,1], expression_data[,1]),]
expression_data_subset_omit_NA = na.omit(expression_data_subset)
