#creating supplementary figures
#May 10th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

#1. Create a table of core genes (252 genes) with their FC in 12hr
core_up_genes = read.table("finding.core.genes/core_genes_for_heatmap/core_upgenes_for_heatmap_function_classified_050817.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data
expression_data_subset = expression_data[match(core_up_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(1,2,9,10,15,16,21,22,39,40,33,34,27,28,45,46,47,48,49,50,51,52)] #for FC - 12hr
write.table(expression_data_subset_with_name, file = "finding.core.genes/core_genes_for_heatmap/supp_core_up_with_FC.txt", quote=F, row.names = F, col.names = T)


#2. Create a table of core genes kinetics (252 genes)
get_Kinetics_with_directions = function(core_genes_address){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    
    heatmap_number = matrix(NA, ncol = (dim(path_data_NS_removed_core_dir_only)[2]-2), nrow = dim(path_data_NS_removed_core_dir_only)[1])
    for (i in 1:dim(path_data_NS_removed_core_dir_only)[1]){ #gene
        for (j in 1:(dim(path_data_NS_removed_core_dir_only)[2]-2)){ #condition
            pattern = as.character(path_data_NS_removed_core_dir_only[i,(j+2)])
            if (pattern == "EE-EE-EE"){ #mark the gene/condition that did not significantly change
                #heatmap_number[i,j] = NA
            }
            
            else{
                #cat("i: ", i, " j: ", j, " j*3: ", j*3, " (j*3)+1: ",(j*3)+1, "\n")
                twelve = expression_data_fc_only[i,j*3]
                thirtysix = expression_data_fc_only[i,(j*3)+1]
                onethirtytwo = expression_data_fc_only[i,(j*3)+2]
                
                ratio=0; onethirtytwo_added =0
                if (twelve > 0 & thirtysix > 0 & onethirtytwo > 0){
                    if (onethirtytwo > twelve & onethirtytwo > thirtysix){
                        pickmax = max(twelve, thirtysix)
                        ratio = round( -abs(onethirtytwo - pickmax)/abs(pickmax), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else{
                        pickmax = max(twelve, thirtysix)
                        ratio = round( abs(onethirtytwo-pickmax)/abs(pickmax), 3)
                        heatmap_number[i,j] = ratio
                    }
                }
                else if (twelve > 0 & thirtysix > 0 & onethirtytwo < 0){
                    pickmax = max(twelve, thirtysix)
                    onethirtytwo_added = pickmax + abs(onethirtytwo)
                    ratio = round(( onethirtytwo_added/pickmax ), 3)
                    heatmap_number[i,j] = ratio
                }
                else if (twelve < 0 & thirtysix > 0 & onethirtytwo > 0){
                    #heatmap_number[i,j] = NA
                    #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                    #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                    #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else {
                    #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                    #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                    #                             heatmap_number[i,j] = ratio
                    #                     }
                } #NA as of 5/1/17
                else if (twelve < 0 & thirtysix < 0 & onethirtytwo > 0){
                    pickmax = max( abs(twelve), abs(thirtysix) )
                    onethirtytwo_added = abs(pickmax) + onethirtytwo
                    ratio = round(( onethirtytwo_added/pickmax ), 3)
                    heatmap_number[i,j] = ratio
                }
                else if (twelve > 0 & thirtysix < 0 & onethirtytwo > 0){
                    #heatmap_number[i,j] = NA
                    #                     if ( abs(twelve) > abs(thirtysix)){
                    #                         ratio = round(( onethirtytwo/twelve ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else{
                    #                         onethirtytwo_added = onethirtytwo + abs(thirtysix)
                    #                         ratio = round(( onethirtytwo_added/ abs(thirtysix) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                      }
                } #NA as of 5/1/17
                else if (twelve > 0 & thirtysix < 0 & onethirtytwo < 0){
                    #heatmap_number[i,j] = NA
                    #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                    #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                    #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else {
                    #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                    #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                } #NA as of 5/1/17
                else if (twelve < 0 & thirtysix > 0 & onethirtytwo < 0){
                    #heatmap_number[i,j] = NA
                    #                     if ( abs(twelve) > abs(thirtysix)){
                    #                         ratio = round(( abs(onethirtytwo)/ abs(twelve) ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                    #                     else{
                    #                         onethirtytwo_added = abs(onethirtytwo) + thirtysix
                    #                         ratio = round(( onethirtytwo_added / thirtysix ), 3)
                    #                         heatmap_number[i,j] = ratio
                    #                     }
                } #NA as of 5/1/17
                else if (twelve < 0 & thirtysix < 0 & onethirtytwo < 0){
                    if ( abs(onethirtytwo) > abs(twelve) & abs(onethirtytwo) > abs(thirtysix)){
                        pickmax = max( abs(twelve), abs(thirtysix) )
                        ratio = round( -abs(onethirtytwo - (-pickmax))/abs(pickmax), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else{
                        pickmax = max( abs(twelve), abs(thirtysix) )
                        ratio = round(( abs(onethirtytwo-(-pickmax)) / pickmax ), 3)
                        heatmap_number[i,j] = ratio
                    }
                }
                
            }
        }
    }
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    
    heatmap_number_with_name = heatmap_number_with_name[match(core_regulated$gene_id, heatmap_number_with_name$gene_id),]
    return(heatmap_number_with_name)
}

#Core upregulated genes:
core_upregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt") #all
core_upregulated_kinetics = core_upregulated_kinetics[,c(1,2,3,11,10,4,5,6,9,8,7)] #change the order to: sw, heatkilleds, ecc15, rest
write.table(core_upregulated_kinetics, file = "clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/supplementary_table_core_upregulated_kinetics.txt", quote=F, row.names = F, col.names = T)

#Core upregulated genes:
core_downregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt") #all
core_downregulated_kinetics = core_downregulated_kinetics[,c(1,2,3,11,10,4,5,6,9,8,7)] #change the order to: sw, heatkilleds, ecc15, rest
write.table(core_downregulated_kinetics, file = "clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/supplementary_table_core_downregulated_kinetics.txt", quote=F, row.names = F, col.names = T)