#Find the recovery genes based on edgeR expression path assignment and sort them by function
#April 25th, 2017
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

library(ggplot2); library(reshape2); library(klaR); library(RColorBrewer)
path_data_NS_removed = read.table("recovery_genes/recovery_pattern_of_genes.txt", header=T)

check_recovery_status_by_category = function(address_to_test_category, include) {
    
    if(include == "Yes"){
        input_data = path_data #includes genes that did not change across all conditions
        colnames(input_data) = colnames(path_data_NS_removed)
    }
    else if(include == "No"){
        input_data= path_data_NS_removed #excludes genes that did not change across all conditions
    }
    
    #this is for function
    #genes_to_test = read.table(address_to_test_category, header=F, sep="\t")
    #genes_to_test[,1] = sapply(genes_to_test[,1], function(x) {sub("FB:","",x)} )
    #genes_to_test = unique(genes_to_test)
    
    #this is for core genes
    genes_to_test = read.table(address_to_test_category, header=T)
    
    input = input_data #excludes genes that did not change across all conditions
    possible_patterns_up = c("Up-Up-EE","Up-EE-EE","EE-Up-EE")
    possible_patterns_down = c("EE-Down-EE","Down-EE-EE", "Down-Down-EE")
    non_recovery_patterns_up = c("EE-EE-Up", "Up-Up-Up", "EE-Up-Up", "Up-EE-Up")
    no_change = c("EE-EE-EE")
    
    genes_to_test_data = c()
    for (i in 1:dim(genes_to_test)[1]){
        genes_to_test_id = as.character(genes_to_test[i,1])
        genes_to_test_name = as.character(genes_to_test[i,2])
        tested_data = input_data[which(input_data$gene.id == genes_to_test_id | input_data$gene.name == genes_to_test_name),]
        if (dim(tested_data)[1] >0){
            genes_to_test_data = rbind(genes_to_test_data, tested_data)
        }
    }
    
    genes_to_test_data_checked = matrix(NA, nrow=dim(genes_to_test_data)[1], ncol=dim(genes_to_test_data)[2])
    for (j in 1:dim(genes_to_test_data)[1]){
        for (k in 3:dim(genes_to_test_data)[2]){
            if (as.character(genes_to_test_data[j,k]) %in% possible_patterns_up){ #if recovered/returned + upregulated
                genes_to_test_data_checked[j,k] = as.character("Yes-up")
            }
            else if (as.character(genes_to_test_data[j,k]) %in% possible_patterns_down){
                genes_to_test_data_checked[j,k] = as.character("Yes-down")
            }
            else if (as.character(genes_to_test_data[j,k]) == no_change) { #if no change
                genes_to_test_data_checked[j,k] = as.character("No_change")
            }
            else if (as.character(genes_to_test_data[j,k]) %in% non_recovery_patterns_up) {
                genes_to_test_data_checked[j,k] = as.character("No-up")
            }
            else {
                genes_to_test_data_checked[j,k] = as.character("No-down")
            }
        }
    }

    data_converted = cbind(genes_to_test_data[,1:2], genes_to_test_data_checked[,3:dim(genes_to_test_data_checked)[2]])
    colnames(data_converted) = c("gene_id", "gene_name", "Sterile Wound","M.luteus", "E.coli", "S.marcescens Type", "E.faecalis", "P.rettgeri", "Ecc15", "E.faecalis heatkilled", "P.rettgeri heatkilled")
    data_converted = data_converted[,c(1:2,4:6,9,8,7,11,10,3)]  #row x 11: changed the order so that its listed based on the virulence level
    
    return(data_converted)
}
cluster_the_genes_by_pattern = function(response_data, sort_by_gene, number_of_cluster){
    if (sort_by_gene == "Yes"){ #cluster the genes (by row)
        cluster.results <-kmodes(response_data[,2:dim(response_data)[2]], number_of_cluster, iter.max = 10, weighted = FALSE)
    }
    else if (sort_by_gene == "No") { #cluster the conditions (by conditions)
        response_data_t = t(response_data)
        cluster.results <-kmodes(response_data_t[2:dim(response_data_t)[1],], number_of_cluster, iter.max = 10, weighted = FALSE)
    }
    
    cluster_info = unlist(cluster.results[1])
    return(cluster_info)
}
plot_the_recovery_status = function(recovery_status_table, sort_by_gene, clustering_info, condition_name, condition_order){
    data = recovery_status_table
    
    if (sort_by_gene == "Yes"){
        data = cbind(data, clustering_info)
        data = data[order(data$clustering_info, decreasing =F),] #sort by clustering info
        data = data[,c(1:(dim(data)[2]-1))] #remove the clustering info
        gene_order = data[,2]
        data$gene_id = droplevels(data$gene_id); data$gene_name = droplevels(data$gene_name)
        data_rearranged <- melt(data, id = c('gene_id', 'gene_name'))
        data_rearranged$variable = factor(data_rearranged$variable, levels= condition_order)
        colnames(data_rearranged) = c("gene_id", "gene_name", "condition", "recovery_pattern")
        data_rearranged$gene_name = factor(data_rearranged$gene_name, levels = gene_order)
    }
    else{ #sort by conditions
        condition_order_with_clustering_info = cbind(condition_order, clustering_info)
        condition_order_with_clustering_info = condition_order_with_clustering_info[order(condition_order_with_clustering_info[,2], decreasing = T),]
        condition_order = condition_order_with_clustering_info[,1]
        
        data$gene_id = droplevels(data$gene_id); data$gene_name = droplevels(data$gene_name)
        data_rearranged <- melt(data, id = c('gene_id', 'gene_name'))
        data_rearranged$variable = factor(data_rearranged$variable, levels= condition_order)
        colnames(data_rearranged) = c("gene_id", "gene_name", "condition", "recovery_pattern")
    }
    
    cols <- c("8" = "red", "4" = "blue", "6" = "darkgreen", "10" = "orange")
    set_colors = c("No_change" = "lightgrey", "No-down" = "#a6611a", "No-up" = "#018571", "Yes-down" = "#dfc27d", "Yes-up" = "#80cdc1")
    plot = ggplot(data_rearranged, aes(condition, gene_name, recovery_pattern)) + geom_tile(aes(fill = recovery_pattern), colour = "white") + scale_fill_manual(values=set_colors) + ggtitle(paste("Recovery pattern for ", condition_name, sep="") ) + theme(text=element_text(size=13), axis.text.x = element_text(angle=90, hjust=1)) + theme(legend.text = element_text(size = 15)) #No_change, No-down, No-up, Yes-down, Yes-up
    return(plot)
}

#if include=yes, it includes genes that did not change across all conditions
#if sort_by_gene=yes, it clusters the genes (by row)
visualize_patterns = function(address_to_test_category, include, sort_by_gene, condition_name, condition_order){
    response_data = check_recovery_status_by_category(address_to_test_category, include)
    number_of_cluster = ceiling(sqrt(dim(response_data)[1])/2)
    clustering_info = cluster_the_genes_by_pattern(response_data, sort_by_gene, number_of_cluster)
    plot_the_recovery_status(response_data, sort_by_gene, clustering_info, condition_name)
} #all put together
condition_order = c("M.luteus","E.coli","S.marcescens Type", "Ecc15", "P.rettgeri", "E.faecalis", "P.rettgeri heatkilled", "E.faecalis heatkilled","Sterile Wound") #ordered by the level of virulence



#Cluster the patterns by genes

#1. Response to Wounding -- 38 genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_wounding_genes_based_on_AmiGO.txt", "No", "Yes","response to wounding", condition_order)


#2. Cuticle development -- 32 genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/cuticle_development_genes_based_on_AmiGO.txt", "No", "Yes", "cuticle development", condition_order)


#3. Response to oxidative stress -- 24 genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_oxidative_stress_genes_based_on_AmiGO.txt", "No", "Yes", "oxidative stress", condition_order)


#4. Humoral immune response/AMP production(lower-level category)
#Humoral immune response -- 53 genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/humoral_immune_response_genes_based_on_AmiGO.txt", "No", "Yes","humoral immune response", condition_order)


#5. Metal ion transport and homeostasis -- 34 genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/metal_ion_merged_genes_based_on_AmiGO.txt", "No", "Yes", "metal ion transport and homeostasis", condition_order)


#6.Toll, JAK-STAT, JNK
#Toll
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/Toll_signaling_pathway_genes_based_on_AmiGO.txt", "No", "Yes", "Toll pathway", condition_order)

#JAK-STAT
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/JAK-STAT_cascade_genes_based_on_AmiGO.txt", "No", "Yes", "JAK-STAT cascade", condition_order)

#JNK
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/JNK_cascade_genes_based_on_AmiGO.txt", "No", "Yes", "JNK cascade", condition_order)


#7. Core genes
#Change the 'check_recovery_status_by_category' function
#upregulated
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", "No", "Yes", "upregulated core genes", condition_order)

#downregulated
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", "No", "Yes", "downregulated core genes", condition_order)



#Look for genes that have different pattern between Ecc15 and E.coli/P.rettgeri (updated on 4/25/2017)

input = path_data_NS_removed #excludes genes that did not change across all conditions
possible_patterns_up = c("Up-Up-EE","Up-EE-EE","EE-Up-EE")
possible_patterns_down = c("EE-Down-EE","Down-EE-EE", "Down-Down-EE")
non_recovery_patterns_up = c("EE-EE-Up", "Up-Up-Up", "EE-Up-Up", "Up-EE-Up")
no_change = c("EE-EE-EE")

genes_to_test_data_checked = matrix(NA, nrow=dim(input)[1], ncol=dim(input)[2])
for (j in 1:dim(input)[1]){
    for (k in 3:dim(input)[2]){
        if (as.character(input[j,k]) %in% possible_patterns_up){ #if recovered/returned + upregulated
            genes_to_test_data_checked[j,k] = as.character("Yes-up")
        }
        else if (as.character(input[j,k]) %in% possible_patterns_down){
            genes_to_test_data_checked[j,k] = as.character("Yes-down")
        }
        else if (as.character(input[j,k]) == no_change) { #if no change
            genes_to_test_data_checked[j,k] = as.character("No_change")
        }
        else if (as.character(input[j,k]) %in% non_recovery_patterns_up) {
            genes_to_test_data_checked[j,k] = as.character("No-up")
        }
        else {
            genes_to_test_data_checked[j,k] = as.character("No-down")
        }
    }
}
data_converted = cbind(input[,1:2], genes_to_test_data_checked[,3:dim(genes_to_test_data_checked)[2]])
colnames(data_converted) = c("gene_id", "gene_name", "Clean.prick", "M.luteus", "E.coli", "S.marcescens", "E.faecalis", "P.rettgeri", "Ecc15","E.faecalis.hk","P.rettgeri.hk")
path_data_NS_removed_converted = data_converted
path_data_NS_removed_chosen_converted = path_data_NS_removed_converted[,c(1:2,9,5,8)]  #Ecc15, E.coli, P.rettgeri
    
#Look for genes that have different pattern between Ecc15 and E.coli/P.rettgeri
path_data_NS_removed_chosen_converted_filtered = c()
for (i in 1:dim(path_data_NS_removed_chosen_converted)[1]){
    #if the pattern in E.coli and P.rett are the same, and the pattern in E.coli and Ecc15 are different
    if (path_data_NS_removed_chosen_converted[i,4] == path_data_NS_removed_chosen_converted[i,5]){
        if (path_data_NS_removed_chosen_converted[i,4] != path_data_NS_removed_chosen_converted[i,3]){
            path_data_NS_removed_chosen_converted_filtered = rbind(path_data_NS_removed_chosen_converted_filtered, path_data_NS_removed_chosen_converted[i,])
        }
    }
}
dim(path_data_NS_removed_chosen_converted_filtered) #216
path_data_NS_removed_chosen_converted_filtered$status = rowSums(path_data_NS_removed_chosen_converted_filtered == "No_change")
path_data_NS_removed_chosen_converted_filtered = path_data_NS_removed_chosen_converted_filtered[which(path_data_NS_removed_chosen_converted_filtered$status ==0),]
dim(path_data_NS_removed_chosen_converted_filtered) #45
write.table(path_data_NS_removed_chosen_converted_filtered, "recovery_genes/recovery_by_cleared_vs_chronic/clustering_info_by_program_manually_corrected_change_only.txt", quote=F, col.names = T, row.names = F)

adjusted_data = read.table("recovery_genes/recovery_by_cleared_vs_chronic/clustering_info_by_program_change_only_man.txt", header=T)
condition_order = c("Ecc15", "E.coli", "P.rettgeri"); gene_order = adjusted_data[,2]
adjusted_data = adjusted_data[,c(1:5)]
adjusted_data$gene_id = droplevels(adjusted_data$gene_id); adjusted_data$gene_name = droplevels(adjusted_data$gene_name)
data_rearranged <- melt(adjusted_data, id = c('gene_id', 'gene_name'))
data_rearranged$variable = factor(data_rearranged$variable, levels= condition_order)
colnames(data_rearranged) = c("gene_id", "gene_name", "condition", "recovery_pattern")
data_rearranged$gene_name = factor(data_rearranged$gene_name, levels = gene_order)
ggplot(data_rearranged, aes(condition, gene_name, recovery_pattern)) + geom_tile(aes(fill = recovery_pattern), colour = "white") + scale_fill_manual(values=c("#a6611a", "#018571", "#dfc27d","#80cdc1")) + ggtitle(paste("Recovery pattern for Ecc15 (cleared) vs E.coli/P.rettgeri (chronic)", sep="") ) + theme(text=element_text(size=12))


#c("lightgrey", "#a6611a", "#018571", "#dfc27d","#80cdc1") = no change, no-down, no-up, yes-down, yes-up



#Kinetics (4/26/2017)
rm(list=ls(all=TRUE))

#Only look at the genes that changed in at least one time point in a given condition (excluding EE-EE-EE).
#Only look at the genes that are consistently upregulated (in the case of core upregulated) or downregulated (in the case of core downregulated). For instance, a core upregulated genes in 7 conditions that is downregulated in another condition has been excluded.

get_Kinetics = function(core_genes_address){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    #expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T) #count
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    #expression_data_count_only = expression_data[,c(1:2,4:24,29:34)] #for count
    
    heatmap_number = matrix(NA, ncol = 9, nrow = dim(path_data_NS_removed_core_dir_only)[1])
    late_start = c()
    for (i in 1:dim(path_data_NS_removed_core_dir_only)[1]){ #gene
        for (j in 1:(dim(path_data_NS_removed_core_dir_only)[2]-2)){ #condition
            pattern = as.character(path_data_NS_removed_core_dir_only[i,(j+2)])
            late_start_row = c()
            #print(pattern)
            
            if (pattern == "EE-EE-EE"){ #no change
                heatmap_number[i,j] = NA
            }
            
            #Mask genes that were core upregulated but then was downregulated in some other conditions
            else if (pattern == "Down-EE-Up" | pattern == "Up-EE-Down"){
                heatmap_number[i,j] = NA
            }
            
            else {
                
                twelve = expression_data_fc_only[i,j*3]
                thirtysix = expression_data_fc_only[i,(j*3)+1]
                onethirtytwo = expression_data_fc_only[i,(j*3)+2]
                
                #if 132hr is maximum: either continuously (Up-Up-Up) or exponentially (EE-EE/Up-Up) increasing or Up-EE-Up
                if ( max(abs(twelve), abs(thirtysix), abs(onethirtytwo)) == abs(onethirtytwo)) { 
                    heatmap_number[i,j] = "late" #NA
                    #late_start_row = matrix(c(as.character(path_data_NS_removed_core_dir_only[i,2]),as.character(colnames(path_data_NS_removed_core_dir_only)[j+2])), 1,2)
                    #late_start = rbind(late_start, late_start_row)
                }
                else {
                    pickmax = max(twelve, thirtysix)
                    ratio = (onethirtytwo/pickmax)*100
                    heatmap_number[i,j] = ratio
                }
            }
        }
    }
    
    #cat("The ones excluded because the expression level at 132h was at its highest: ", late_start)
    
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    return(heatmap_number_with_name)
}

core_upregulated_kinetics = get_Kinetics("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt")
core_downregulated_kinetics = get_Kinetics("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt")











