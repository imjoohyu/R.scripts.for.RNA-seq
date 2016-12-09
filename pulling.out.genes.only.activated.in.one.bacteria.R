#EdgeR Commands for Mega RNA-seq for pulling out genes that only [insert the bacteria name] activated
#Date: October 8th, 2015 (Updated on Feb 5th, 2016, further updated on 11/8/2016, 12/7-9/2016)
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data_table <- read.table("edgeR_basic_comparison_pval_at_least_one_sig.txt",header=T)

#Prep: Assign the direction by converting the FC information to Up, Down, or Same depending on their degree of significance
length.of.table <- as.numeric(dim(data_table)[1])
width.of.table <- as.numeric(dim(data_table)[2])
indicator <- NULL

for (k in 1:length.of.table){ #1, 2, 3, ... 2589
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        indicator <- isTRUE(data_table[k,m] > 0.05) #cutoff: FDR of 5%
        if (indicator == TRUE) { #If the case is NOT significant,
            data_table[k,m] = "N"
        }
        else {
            data_table[k,m] = "Y"
        }
    }
}
for (i in 1:length.of.table){ #1, 2, 3, ... 2589
    for (s in seq(3, width.of.table, 2)){ #3, 5, 7, ... 63
        
        if (data_table[i,s+1] == "Y"){ #if significant
            indicator <- isTRUE(data_table[i,s] > 0) #indicator shows that the direction is positive
            #cat("indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the case is Up-DEG
                data_table[i,s] = "Up"
            }
            else { #If the caseis Down-DEG
                data_table[i,s] = "Down"
            }
        }
        else { #if not significant
            data_table[i,s] = "EE"
        }
        
    }
}


####
#1. From the data.table, find the genes that were significant in a given condition
#2. Now with the list of genes that were significant in a given condition, only select genes that have no significant expression in any other infection conditions.
condition_list = c(3,9,15,21,27,53,33,59,39,45,47,49,51) #12hr time point
condition_name_list = c("SterileWound","M.luteus","E.coli","S.marcescens","E.faecalis","E.faecalis.hk","P.rettgeri","P.rettgeri.hk","Ecc15","S.aureus","P.sneebia","S.marcescens_Db11","P.entomophila")
percentage_of_unique_genes=matrix(NA, ncol=4, nrow=length(condition_list))

pull_out_genes_only_activated_in_one_condition = function(direction){
    for (i in 1:length(condition_list)){
        condition_number = condition_list[i]
        subset_data_table = subset(data_table, grepl(direction,data_table[,condition_number])) #for a given condition, pick the genes that had the chosen direction.
        total_number_of_DEGs_in_this_condition = dim(subset_data_table)[1]
        condition_list_other_than_chosen_condition = condition_list[!condition_list == condition_number]
        
        cat('\n')
        cat("condition_list_other_than_chosen_condition: ", condition_list_other_than_chosen_condition, '\n')
        subset_data_table_rest = subset_data_table[,c(1,2,condition_list_other_than_chosen_condition)]
        
        subset_data_table_cross_check_with_other_conditions=c()
        for (j in 1:dim(subset_data_table_rest)[1]){
            gene = subset_data_table_rest[j,c(3:14)]
            count = sum(gene != direction)  #if the number of 'EE' or the opposite of direction is = 12 -> this gene is unique
            if (count == 12) {#only when this gene is unique
                subset_data_table_cross_check_with_other_conditions = rbind(subset_data_table_cross_check_with_other_conditions, subset_data_table[j,c(1,2,condition_list)])
            }
        }
        
        subset_data_table_cross_check_with_other_conditions_final = subset_data_table_cross_check_with_other_conditions
        if (is.null(subset_data_table_cross_check_with_other_conditions) == TRUE){
            total_number_of_unique_DEGs_in_this_condition = 0
        }
        else{
            total_number_of_unique_DEGs_in_this_condition = dim(subset_data_table_cross_check_with_other_conditions_final)[1]
        }
        
        cat("Infection condition: ", condition_name_list[i], '\n')
        cat("total_number_of_DEGs_in_this_condition: ", total_number_of_DEGs_in_this_condition, '\n')
        cat("total_number_of_unique_DEGs_in_this_condition: ", total_number_of_unique_DEGs_in_this_condition, '\n')
        cat("percentage of unique genes out of total DEGs: ", as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition, '\n'))
        
        percentage_of_unique_genes[i,1] = condition_name_list[i]
        percentage_of_unique_genes[i,2] = total_number_of_unique_DEGs_in_this_condition
        percentage_of_unique_genes[i,3] = total_number_of_DEGs_in_this_condition
        percentage_of_unique_genes[i,4] = as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition)
        
    }
    
    colnames(percentage_of_unique_genes) = c("condition","total_number_of_unique_DEGs_in_this_condition", "total_number_of_DEGs_in_this_condition", "percentage_of_unique_genes")
    if (direction == "Up"){
        write.table(percentage_of_unique_genes, file = "../specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_percentage_of_unique_upreg_genes_including_cleanprick_and_heatkilled.txt", col.names = T, row.names = F, quote=F)
    }
    if (direction == "Down"){
        write.table(percentage_of_unique_genes, file = "../specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_percentage_of_unique_downreg_genes_including_cleanprick_and_heatkilled.txt", col.names = T, row.names = F, quote=F)
    }
    
    return(percentage_of_unique_genes)
}
up_genes = pull_out_genes_only_activated_in_one_condition("Up")
down_genes = pull_out_genes_only_activated_in_one_condition("Down")


#3. Plot this percentage in a bar chart
par(mfrow = c(2,1))
Draw_a_graph = function(data, direction){
    color_list = c("azure4","lightpink", "darkseagreen2","lightseagreen", "coral1","dark salmon","dodgerblue2","steelblue1","lightskyblue","orangered3", "blue", "navy", "purple3")
    
    percentage_table = c(data[,4]); percentage_list = as.numeric(percentage_table); percentage_list = percentage_list*100
    
    if (direction == "Up"){
        barplot(percentage_list, main="Percentage of unique-DEGs per infection condition", ylab = "Percentage of unique Up-DEGs (%)", names.arg=condition_name_list, col= color_list, ylim=c(0,50), las=2)
        text(0.72, 5,"12.1%"); text(1.88, 5, "15.4%"); text(3.1, 5, "13.5%"); text(4.3, 5, "8.9%")
        text(5.5, 5, "20.7%"); text(6.7, 5, "10.6%"); text(7.9, 5, "19.6%"); text(9.1, 5, "17.7%")
        text(10.3, 5, "9.9%"); text(11.5, 5, "37.4%"); text(12.7, 5, "10.8%"); text(13.9, 5, "15.5%",col="white")
        text(15.1, 5, "40.2%")
    }
    else{
        barplot(percentage_list, ylab = "Percentage of unique Down-DEGs (%)", names.arg=condition_name_list, col= color_list, ylim=c(50,0), las=2)
        text(0.72, 5,"0%"); text(1.88, 5, "16.0%"); text(3.1, 5, "8.6%"); text(4.3, 5, "6.4%")
        text(5.5, 5, "10.1%"); text(6.7, 5, "20.0%"); text(7.9, 5, "20.2%"); text(9.1, 5, "39.1%")
        text(10.3, 5, "11.9%"); text(11.5, 5, "47.7%"); text(12.7, 5, "20.6%"); text(13.9, 5, "19.8%",col="white")
        text(15.1, 5, "29.4%")
    }
}
Draw_a_graph(up_genes, "Up"); Draw_a_graph(down_genes, "Down")


#4. Create a dendrogram based on the original data
data_table <- read.table("edgeR_basic_comparison_pval_at_least_one_sig.txt",header=T)
library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("red","black","green"))(n = 299)

#Dendrogram based on every sample
odd_index = seq(3,64,2)
data_table_FC_only = data_table[,c(odd_index)]; rownames(data_table_FC_only) = data_table[,1]
heatmap.2(as.matrix(data_table_FC_only), Rowv=FALSE, density.info="none", dendrogram="column", trace="none", scale=c("column"), col=my_palette)

#Dendrogram based on 12hr sample >> This pattern goes well with the PCA pattern as expected (12/9/2016)
data_table_FC_only_12hr = data_table_FC_only[,c(1,4,7,10,13,16,19,22,23,24,25,26,29)] #12hr sample only
heatmap.2(as.matrix(data_table_FC_only_12hr), Rowv=FALSE, density.info="none", dendrogram="column", trace="none", scale=c("column"), col=my_palette, srtCol=45)


#5. Change the order of the data table and plot a bar chart in the order of the dendrogram
up_genes_rearranged = rbind(up_genes[10,],up_genes[7,],up_genes[13,],up_genes[1,],up_genes[9,],up_genes[6,],up_genes[8,],up_genes[5,],up_genes[2,],up_genes[12,],up_genes[11,],up_genes[4,],up_genes[3,])
down_genes_rearranged = rbind(down_genes[10,],down_genes[7,],down_genes[13,],down_genes[1,],down_genes[9,],down_genes[6,],down_genes[8,],down_genes[5,],down_genes[2,],down_genes[12,],down_genes[11,],down_genes[4,],down_genes[3,])


par(mfrow = c(2,1))
Draw_a_graph_rev = function(data, direction){
    color_list_2 = c("orangered3","dodgerblue2","purple3","azure4","lightskyblue","dark salmon","steelblue1","coral1","lightpink","navy","blue","lightseagreen","darkseagreen2")
    
    percentage_table = c(data[,4]); percentage_list = as.numeric(percentage_table); percentage_list = percentage_list*100
    coor=c(0.72, 1.88, 3.1, 4.3, 5.5, 6.7, 7.9, 9.1, 10.3, 11.5, 12.7, 13.9, 15.1)
    
    if (direction == "Up"){
        barplot(percentage_list, main="Percentage of unique-DEGs per infection condition", ylab = "Percentage of unique Up-DEGs (%)", names.arg=data[,1], col= color_list_2, ylim=c(0,50), las=2)
        
        for (m in 1:length(percentage_list)){
            if (m == 10){
                text(coor[m], 5, paste(round(percentage_list[m],1),"%",sep=""), col="white")
            }
            else {
                text(coor[m], 5, paste(round(percentage_list[m],1),"%",sep=""))
            }
        }
    }
    else{
        barplot(percentage_list, ylab = "Percentage of unique Down-DEGs (%)", names.arg=data[,1], col= color_list_2, ylim=c(50,0), las=2)
        for (m in 1:length(percentage_list)){
            if (m == 10){
                text(coor[m], 5, paste(round(percentage_list[m],1),"%",sep=""), col="white")
            }
            else {
                text(coor[m], 5, paste(round(percentage_list[m],1),"%",sep=""))
            }
        }
    }
}
Draw_a_graph_rev(up_genes_rearranged, "Up"); Draw_a_graph_rev(down_genes_rearranged, "Down")






