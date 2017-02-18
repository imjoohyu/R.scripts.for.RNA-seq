#EdgeR Commands for Mega RNA-seq for pulling out genes that only uniquely activated by [insert the bacteria name]
#Date: October 8th, 2015 (Updated on Feb 5th, 2016, further updated on 11/8/2016, 12/7-9/2016)
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data_table = read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)

#Prep: Assign the direction by converting the FC information to Up, Down, or Same depending on their degree of significance
length.of.table <- as.numeric(dim(data_table)[1])
width.of.table <- as.numeric(dim(data_table)[2])
indicator <- NULL

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
#Definition of an uniquely expressed gene:
#Method I. a gene that changes expression significantly in one and only one infection condition at 12hrs post infection.
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
        
        write.table(subset_data_table_cross_check_with_other_conditions_final, file=paste("../specific.comparisons/genes_uniquely_regulated_by_each_condition/uniquely_",direction,"regulated_by_",condition_name_list[i],".txt",sep=""), col.names=T, row.names=F, quote=F)
        
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



#Method II. a gene that changes expression significantly in one and only one infection condition regardless of time points (collapsing time).
#1. From the data.table, find the genes that were significant in a given condition
#2. Now with the list of genes that were significant in a given condition, only select genes that have no significant expression in any other infection conditions.
condition_list = c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,53,55,57,33,35,37,59,61,63,39,41,43,45,47,49,51) #all time points
condition_name_list = c("SterileWound","M.luteus","E.coli","S.marcescens","E.faecalis","E.faecalis.hk","P.rettgeri","P.rettgeri.hk","Ecc15","S.aureus","P.sneebia","S.marcescens_Db11","P.entomophila")
percentage_of_unique_genes=matrix(NA, ncol=4, nrow=length(condition_name_list))

pull_out_genes_only_activated_in_one_condition_ignoring_time = function(direction){
    for (k in 1:length(condition_name_list)){
        condition_number = NULL; condition_list_other_than_chosen_condition = NULL
        if (k < 10){ #for conditions with all three time points
            condition_number = condition_list[(k*3-2):(k*3)]
            #for a given condition, pick the genes that had the chosen direction regardless of time points (pick it if it is significant in at least one time point)
            subset_data_table = data_table[with(data_table, grepl(direction,data_table[,condition_number[1]]) | grepl(direction,data_table[,condition_number[2]]) | grepl(direction,data_table[,condition_number[3]])), ] 
            total_number_of_DEGs_in_this_condition = dim(subset_data_table)[1]
            condition_list_other_than_chosen_condition = condition_list[which(!condition_list == unique(condition_number, condition_list))] #generates warning messages
            
        }
        else{ #conditions with 1 time point
            condition_number = condition_list[(k+18)]
            subset_data_table = subset(data_table, grepl(direction,data_table[,condition_number]))
            total_number_of_DEGs_in_this_condition = dim(subset_data_table)[1]
            condition_list_other_than_chosen_condition = condition_list[which(!condition_list == unique(condition_number,condition_list))] #generates warning messages
        }
        
        cat('\n')
        cat("condition_being_checked: ", condition_name_list[k], '\n')
        cat("condition_list: ", condition_list, '\n')dim
        cat("condition_list_other_than_chosen_condition: ", condition_list_other_than_chosen_condition, '\n')
        subset_data_table_rest = subset_data_table[,c(1,2,condition_list_other_than_chosen_condition)]
        cat("names: ", colnames(subset_data_table_rest), '\n')
        
        subset_data_table_cross_check_with_other_conditions=c()
        for (m in 1:dim(subset_data_table_rest)[1]){
            count=NULL
            if (k <10){
                gene = subset_data_table_rest[m, c(3:30)]
            }
            else{
                gene = subset_data_table_rest[m, c(3:32)]
            }
            
            count = sum(gene == direction)  #count the number of direction
            if (count == 0) {#If there's no other condition that has this direction
                subset_data_table_cross_check_with_other_conditions = rbind(subset_data_table_cross_check_with_other_conditions, subset_data_table[m,c(1,2,condition_list)])
            }
        }
        
        subset_data_table_cross_check_with_other_conditions_final = subset_data_table_cross_check_with_other_conditions
        if (is.null(subset_data_table_cross_check_with_other_conditions) == TRUE){
            total_number_of_unique_DEGs_in_this_condition = 0
        }
        else{
            total_number_of_unique_DEGs_in_this_condition = dim(subset_data_table_cross_check_with_other_conditions_final)[1]
        }
        
        cat("Infection condition: ", condition_name_list[k], '\n')
        cat("total_number_of_DEGs_in_this_condition: ", total_number_of_DEGs_in_this_condition, '\n')
        cat("total_number_of_unique_DEGs_in_this_condition: ", total_number_of_unique_DEGs_in_this_condition, '\n')
        cat("percentage of unique genes out of total DEGs: ", as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition, '\n'))
        
        percentage_of_unique_genes[k,1] = condition_name_list[k]
        percentage_of_unique_genes[k,2] = total_number_of_unique_DEGs_in_this_condition
        percentage_of_unique_genes[k,3] = total_number_of_DEGs_in_this_condition
        percentage_of_unique_genes[k,4] = as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition)
        
        write.table(subset_data_table_cross_check_with_other_conditions_final, file=paste("../specific.comparisons/genes_uniquely_regulated_by_each_condition/uniquely_",direction,"regulated_by_",condition_name_list[k],"_ignoring_time.txt",sep=""), col.names=T, row.names=F, quote=F)
        
    }
    
    colnames(percentage_of_unique_genes) = c("condition","total_number_of_unique_DEGs_in_this_condition", "total_number_of_DEGs_in_this_condition", "percentage_of_unique_genes")
    if (direction == "Up"){
        write.table(percentage_of_unique_genes, file = "../specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_percentage_of_unique_upreg_genes_including_cleanprick_and_heatkilled_ignoring_time.txt", col.names = T, row.names = F, quote=F)
    }
    if (direction == "Down"){
        write.table(percentage_of_unique_genes, file = "../specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_percentage_of_unique_downreg_genes_including_cleanprick_and_heatkilled_ignoring_time.txt", col.names = T, row.names = F, quote=F)
    }
    
    return(percentage_of_unique_genes)
}
up_genes = pull_out_genes_only_activated_in_one_condition_ignoring_time("Up")
down_genes = pull_out_genes_only_activated_in_one_condition_ignoring_time("Down")



#Graphics:

# #3. Plot this percentage in a bar chart
# par(mfrow = c(2,1))
# Draw_a_graph = function(data, direction){
#     color_list = c("azure4","lightpink", "darkseagreen2","lightseagreen", "coral1","dark salmon","dodgerblue2","steelblue1","lightskyblue","orangered3", "blue", "navy", "purple3")
#     
#     percentage_table = c(data[,4]); percentage_list = as.numeric(percentage_table); percentage_list = percentage_list*100
#     
#     if (direction == "Up"){
#         barplot(percentage_list, main="Percentage of unique-DEGs per infection condition", ylab = "Percentage of unique Up-DEGs (%)", names.arg=condition_name_list, col= color_list, ylim=c(0,50), las=2)
#         text(0.72, 5,"12.1%"); text(1.88, 5, "15.4%"); text(3.1, 5, "13.5%"); text(4.3, 5, "8.9%")
#         text(5.5, 5, "20.7%"); text(6.7, 5, "10.6%"); text(7.9, 5, "19.6%"); text(9.1, 5, "17.7%")
#         text(10.3, 5, "9.9%"); text(11.5, 5, "37.4%"); text(12.7, 5, "10.8%"); text(13.9, 5, "15.5%",col="white")
#         text(15.1, 5, "40.2%")
#     }
#     else{
#         barplot(percentage_list, ylab = "Percentage of unique Down-DEGs (%)", names.arg=condition_name_list, col= color_list, ylim=c(50,0), las=2)
#         text(0.72, 5,"0%"); text(1.88, 5, "16.0%"); text(3.1, 5, "8.6%"); text(4.3, 5, "6.4%")
#         text(5.5, 5, "10.1%"); text(6.7, 5, "20.0%"); text(7.9, 5, "20.2%"); text(9.1, 5, "39.1%")
#         text(10.3, 5, "11.9%"); text(11.5, 5, "47.7%"); text(12.7, 5, "20.6%"); text(13.9, 5, "19.8%",col="white")
#         text(15.1, 5, "29.4%")
#     }
# }
# Draw_a_graph(up_genes, "Up"); Draw_a_graph(down_genes, "Down")


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
        barplot(percentage_list, main="Percentage of unique-DEGs per infection condition", ylab = "Percentage of unique Up-DEGs (%)", names.arg=data[,1], col= color_list_2, ylim=c(0,100), las=2)
        
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
        barplot(percentage_list, ylab = "Percentage of unique Down-DEGs (%)", names.arg=data[,1], col= color_list_2, ylim=c(100,0), las=2)
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


#Added on 2/16/2017 (Thur)
#6.Plot the number of DEGs and overlay the number of unique DEGs on top.
dim(data_table) #input data with the direction

condition_name_list = c("SterileWound","M.luteus","E.coli","S.marcescens","E.faecalis","E.faecalis.hk","P.rettgeri","P.rettgeri.hk","Ecc15","S.aureus","P.sneebia","S.marcescens_Db11","P.entomophila")
condition_list = c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,53,55,57,33,35,37,59,61,63,39,41,43,45,47,49,51) #all time points

Get_the_number_of_DEGs_and_unique_genes = function(input_data, direction){
    table_of_DEG_and_unique_DEG = c()
    
    for (i in 1:length(condition_name_list)){ #for each condition
        
        unique_genes = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/genes_uniquely_regulated_by_each_condition/uniquely_",direction,"regulated_by_",condition_name_list[i],"_ignoring_time.txt",sep=""), header=T)
        
        unique_genes_only = unique_genes[,1]
        tp12 = NULL; tp36 = NULL; tp55 = NULL; data_to_add = NULL
        
        if (i < 10){ #for conditions with all three time points
            tp12 = input_data[which(input_data[,condition_list[i*3-2]] == direction),]
            tp36 = input_data[which(input_data[,condition_list[i*3-1]] == direction),]
            tp55 = input_data[which(input_data[,condition_list[i*3]] == direction),]
        
            tp12_intersect = intersect(unique_genes_only, tp12[,1])
            tp36_intersect = intersect(unique_genes_only, tp36[,1])
            tp55_intersect = intersect(unique_genes_only, tp55[,1])
            
            data_to_add = c(condition_name_list[i], dim(tp12)[1], dim(tp36)[1], dim(tp55)[1], length(tp12_intersect), length(tp36_intersect), length(tp55_intersect))
        }
        else { #for conditions with one time point only
            tp12 = input_data[which(input_data[,condition_list[i+18]] == direction),]
            tp12_intersect = intersect(unique_genes_only, tp12[,1])
            
            data_to_add = c(condition_name_list[i], dim(tp12)[1], 0, 0, length(tp12_intersect), 0, 0)
        }
        table_of_DEG_and_unique_DEG = rbind(table_of_DEG_and_unique_DEG, data_to_add)
    }
    return(table_of_DEG_and_unique_DEG)
}
upregulated_uniq_DEG = Get_the_number_of_DEGs_and_unique_genes(data_table, "Up"); colnames(upregulated_uniq_DEG) = c("condition","12hr","36hr","132hr","12hr_uniq","36hr_uniq","132hr_uniq")
downregulated_uniq_DEG = Get_the_number_of_DEGs_and_unique_genes(data_table, "Down"); colnames(downregulated_uniq_DEG) = c("condition","12hr","36hr","132hr","12hr_uniq","36hr_uniq","132hr_uniq") #has an error msg when i=1 because there's no unique gene. But the code works.
write.table(upregulated_uniq_DEG, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_DEG_and_unique_DEG_Upregulated.txt", quote=F, row.names = F, col.names=T)
write.table(downregulated_uniq_DEG, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_DEG_and_unique_DEG_Downregulated.txt", quote=F, row.names = F, col.names=T)


format_the_data = function(input_data){
    input_data_df = data.frame(input_data)
    input_data_df[,1] = as.character(input_data_df[,1])
    for (m in 2:7){
        input_data_df[,m] <- as.numeric(as.character(input_data_df[,m]))
    }
    output = c()
    for (j in 1:length(condition_name_list)){
        tp12 = NULL; tp36 = NULL; tp55 = NULL; data_to_add = NULL
        if (j<10){
            tp12 = c(condition_name_list[j],colnames(input_data_df)[2],"DEG", input_data_df[j,2]-input_data_df[j,5])
            tp12_2 = c(condition_name_list[j],colnames(input_data_df)[2],"Unique", input_data_df[j,5])
            tp36 = c(condition_name_list[j],colnames(input_data_df)[3],"DEG", input_data_df[j,3]-input_data_df[j,6])
            tp36_2 = c(condition_name_list[j],colnames(input_data_df)[3],"Unique", input_data_df[j,6])
            tp55 = c(condition_name_list[j],colnames(input_data_df)[4],"DEG", input_data_df[j,4]-input_data_df[j,7])
            tp55_2 = c(condition_name_list[j],colnames(input_data_df)[4],"Unique", input_data_df[j,7])
            output = rbind(output, tp12, tp12_2, tp36, tp36_2, tp55, tp55_2)
        }
        else {
            tp12 = c(condition_name_list[j],colnames(input_data_df)[2],"DEG", input_data_df[j,2]-input_data_df[j,5])
            tp12_2 = c(condition_name_list[j],colnames(input_data_df)[2],"Unique", input_data_df[j,5])
            tp36 = c(condition_name_list[j],colnames(input_data_df)[3],"DEG",0)
            tp36_2 = c(condition_name_list[j],colnames(input_data_df)[3],"Unique",0)
            tp55 = c(condition_name_list[j],colnames(input_data_df)[4],"DEG",0)
            tp55_2 = c(condition_name_list[j],colnames(input_data_df)[4],"Unique",0)
            output = rbind(output, tp12, tp12_2, tp36, tp36_2, tp55, tp55_2)        
        }
    }
    
    return(output)
} #gets an error msg but it works.
upregulated_uniq_DEG_formatted = as.data.frame(format_the_data(upregulated_uniq_DEG)); colnames(upregulated_uniq_DEG_formatted) = c("Condition","Time","Type","Count")
upregulated_uniq_DEG_formatted[,2] = sub("X12hr", "12hr", x = upregulated_uniq_DEG_formatted[,2] )
upregulated_uniq_DEG_formatted[,2] = sub("X36hr", "36hr", x = upregulated_uniq_DEG_formatted[,2] )
upregulated_uniq_DEG_formatted[,2] = sub("X132hr", "132hr", x = upregulated_uniq_DEG_formatted[,2] )
downregulated_uniq_DEG_formatted = as.data.frame(format_the_data(downregulated_uniq_DEG)); colnames(downregulated_uniq_DEG_formatted) =c("Condition","Time","Type","Count")
downregulated_uniq_DEG_formatted[,2] =sub("X12hr", "12hr", x = downregulated_uniq_DEG_formatted[,2] )
downregulated_uniq_DEG_formatted[,2] =sub("X36hr", "36hr", x = downregulated_uniq_DEG_formatted[,2] )
downregulated_uniq_DEG_formatted[,2] =sub("X132hr", "132hr", x = downregulated_uniq_DEG_formatted[,2] )


library(ggplot2); library(grid); library(gridExtra)
grid.newpage(); pushViewport(viewport(layout = grid.layout(2,1 )))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

Draw_a_plot = function(input_data, direction){
    #condition_order = c("SterileWound","M.luteus","E.coli","S.marcescens","E.faecalis","E.faecalis.hk","P.rettgeri","P.rettgeri.hk","Ecc15","S.aureus","P.sneebia","S.marcescens_Db11","P.entomophila")
    condition_order = c("S.aureus","P.rettgeri","P.entomophila","SterileWound","Ecc15","E.faecalis.hk","P.rettgeri.hk","E.faecalis","M.luteus","S.marcescens_Db11","P.sneebia","S.marcescens","E.coli") #order by dendrogram
    time_order = c("12hr", "36hr", "132hr")
    if (direction == "Up"){
        input_data[,4] = as.numeric(as.character(input_data[,4]))
        input_data$Condition = factor(input_data$Condition, levels= condition_order)
        input_data$Time = factor(input_data$Time, levels= time_order)
        plot = ggplot(input_data, aes(x=Time,y=Count,fill=Type)) + geom_bar(stat = "identity",color="white") + facet_wrap(~Condition,nrow=1) + scale_fill_manual(values = c("black", "orange")) + scale_y_continuous(limits=c(0,750)) + theme_bw(base_size=14)
    }
    else{
        input_data[,4] = as.numeric(as.character(input_data[,4]))
        input_data$Condition = factor(input_data$Condition, levels= condition_order)
        input_data$Time = factor(input_data$Time, levels= time_order)
        plot = ggplot(input_data, aes(x=Time,y=Count,fill=Type)) + geom_bar(stat = "identity",color="white") + facet_wrap(~Condition,nrow=1) + scale_fill_manual(values = c("black", "orange")) + scale_y_continuous(limits=c(0,750)) + theme_bw(base_size=14) + scale_y_reverse()
    }
    return(plot)
}
plot1 = Draw_a_plot(upregulated_uniq_DEG_formatted, "Up")
plot2 = Draw_a_plot(downregulated_uniq_DEG_formatted, "Down")

print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(2, 1))
