#What are the genes were upregulated by some bacteria but downregulated by others?
#March 7th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))

#Read in data
original_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)

#Assign the direction
assign_direction = function(data){
    length.of.table <- as.numeric(dim(data)[1])
    width.of.table <- as.numeric(dim(data)[2])
    indicator <- NULL
    
    for (i in 1:length.of.table){
        for (s in seq(3, width.of.table, 2)){
            
            if (data[i,s+1] == "Y"){ #if significant
                indicator <- isTRUE(data[i,s] > 0) #indicator shows that the direction is positive
                if (indicator == TRUE) { #If the case is Up-DEG
                    data[i,s] = "Up"
                }
                else { #If the caseis Down-DEG
                    data[i,s] = "Down"
                }
            }
            else { #if not significant
                data[i,s] = "EE"
            }
            
        }
    }
    return(data)
}
direction_data = assign_direction(original_data)

#Look for genes that are upregulated by some conditions but are downregulated by others
#https://stat.ethz.ch/pipermail/r-help/2005-February/065381.html
#select rows with 'Up's
direction_data_subset_ups = subset(direction_data, apply(direction_data, 1, function(x){any(x == "Up")}))
dim(direction_data_subset_ups)
#select rows with "Down's
direction_data_subset_downs = subset(direction_data, apply(direction_data, 1, function(x){any(x == "Down")}))
dim(direction_data_subset_downs)
#of the ones with 'Up's, select rows with 'Down's
direction_data_subset_ups_or_downs = subset(direction_data_subset_ups, apply(direction_data_subset_ups, 1, function(x){any(x == "Down")}))
dim(direction_data_subset_ups_or_downs)

write.table(direction_data_subset_ups_or_downs, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others.txt", quote=F, row.names = F, col.names = T)

#Draw out the heatmap-like figure
library(ggplot2); library(reshape2); library(stringr)

direction_data_subset_ups_or_downs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others_w_functions.txt", header=T, sep="\t") #function added
direction_data_subset_ups_or_downs_FC_only = direction_data_subset_ups_or_downs[,c(1:3,seq(4,65,2))]
direction_data_subset_ups_or_downs_FC_only = direction_data_subset_ups_or_downs_FC_only[order(direction_data_subset_ups_or_downs_FC_only$gene_function),]

data = direction_data_subset_ups_or_downs_FC_only[,c(1:2,4:34)]
data$gene_id = droplevels(data$gene_id); data$gene_name = droplevels(data$gene_name)
data_rearranged <- melt(data, id = c('gene_id', 'gene_name'))
data_rearranged$variable = data_rearranged$variable %>% str_replace("log2FC.", "") #simplify the condition names
condition_order = unique(data_rearranged$variable)
data_rearranged$variable = factor(data_rearranged$variable, levels= c(condition_order))

table(direction_data_subset_ups_or_downs_FC_only[,3]) # 2, 76, 1, 4
function_order = c(rep(c("immunity"),2), rep(c("metabolism"),76), "neuron-related", rep(c("transmembrane transport"),4),rep(c("NA"),88))
data_rearranged_function_added = cbind(data_rearranged, rep(function_order,31))
colnames(data_rearranged_function_added) = c("gene_id", "gene_name", "condition", "direction", "gene_function")

#from: http://stackoverflow.com/questions/33271977/order-of-columns-and-rows-with-ggplot2-tile
gene_order = c("immunity", "metabolism", "neuron-related", "transmembrane transport", "NA") #set vector of levels I want
data_rearranged_function_added$gene_function =factor(data_rearranged_function_added$gene_function, levels=gene_order)
gene_name_order = unique(data_rearranged$gene_name)
data_rearranged_function_added$gene_name =factor(data_rearranged_function_added$gene_name, levels=gene_name_order)

plot = ggplot(data_rearranged_function_added, aes(condition, gene_name, direction)) + geom_tile(aes(fill = direction), colour = "white") + scale_fill_manual(values=c("#d8b365", "lightgrey", "#5ab4ac")) + ggtitle("Genes upregulated in certain conditions and downregulateed in other conditions ") + theme(text=element_text(size=15)) + ylab(data_rearranged_function_added$gene_name)





#Look for genes that are upregulated by some conditions but are downregulated by others (12hr only)
direction_data_12hr_only = direction_data[,c(1,2,3,9,15,21,27,33,39,45,47,49,51,53,59)]
direction_data_12hr_only_up = subset(direction_data_12hr_only, apply(direction_data_12hr_only, 1, function(x){any(x == "Up")}))
dim(direction_data_12hr_only_up)
direction_data_12hr_only_up_and_down = subset(direction_data_12hr_only_up, apply(direction_data_12hr_only_up, 1, function(x){any(x == "Down")}))
dim(direction_data_12hr_only_up_and_down)
write.table(direction_data_12hr_only_up_and_down, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others_12hr_only.txt", quote=F, row.names = F, col.names = T)

