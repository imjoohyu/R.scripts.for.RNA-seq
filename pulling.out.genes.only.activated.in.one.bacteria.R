#EdgeR Commands for Mega RNA-seq for pulling out genes that only [insert the bacteria name] activated
#Date: October 8th, 2015 (Updated on Feb 5th, 2016, further updated on 11/8/2016)
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data_table <- read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt",header=T)

####
#1. For fair comparisons, I only looked at 12hr time point and live infection
#2. This method counts both upregulated genes and downregulated genes together


#1. From the data.table, find the genes that were significant in a given condition
#2. Now with the list of genes that were significant in a given condition, only select genes that have no significant expression in any other infection conditions.
condition_list = c(10,16,22,28,34,40,46,48,50,52) #12hr time point
condition_name_list = c("M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","S.aureus","P.sneebia","S.marcescens_Db11","P.entomophila")
percentage_of_unique_genes=matrix(NA, ncol=4, nrow=length(condition_list))


for (i in 1:length(condition_list)){
    condition_number = condition_list[i]
    subset_data_table = subset(data_table, grepl("Y",data_table[,condition_number]))
    total_number_of_DEGs_in_this_condition = dim(subset_data_table)[1]
    condition_list_other_than_chosen_condition = condition_list[!condition_list == condition_number]
    
    cat('\n')
    cat("condition_list_other_than_chosen_condition: ", condition_list_other_than_chosen_condition, '\n')
    subset_data_table_rest = subset_data_table[,c(1,2,condition_list_other_than_chosen_condition)]
    
    subset_data_table_cross_check_with_other_conditions=c()
    for (j in 1:dim(subset_data_table_rest)[1]){
        gene = subset_data_table_rest[j,c(3:11)]
        count = sum(gene == "N")  #if the number of 'N' is = 9 -> this gene is unique
        if (count == 9) {#only when this gene is unique
            subset_data_table_cross_check_with_other_conditions = rbind(subset_data_table_cross_check_with_other_conditions, subset_data_table[j,c(1,2,condition_list)])
        }
    }
    
    subset_data_table_cross_check_with_other_conditions_final = subset_data_table_cross_check_with_other_conditions
    total_number_of_unique_DEGs_in_this_condition = dim(subset_data_table_cross_check_with_other_conditions_final)[1]
    cat("Infection condition: ", condition_name_list[i], '\n')
    cat("total_number_of_DEGs_in_this_condition: ", total_number_of_DEGs_in_this_condition, '\n')
    cat("total_number_of_unique_DEGs_in_this_condition: ", total_number_of_unique_DEGs_in_this_condition, '\n')
    cat("percentage of unique genes out of total DEGs: ", as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition, '\n'))
    
    percentage_of_unique_genes[i,1] = condition_name_list[i]
    percentage_of_unique_genes[i,2] = total_number_of_unique_DEGs_in_this_condition
    percentage_of_unique_genes[i,3] = total_number_of_DEGs_in_this_condition
    percentage_of_unique_genes[i,4] = as.numeric(total_number_of_unique_DEGs_in_this_condition/total_number_of_DEGs_in_this_condition)
    
    #write.table(subset_data_table_cross_check_with_other_conditions_final, file=paste("../specific.comparisons/genes_uniquely_regulated_by_each_condition/uniquely_regulated_by_", condition_name_list[i],".txt", sep=""), row.names = F, col.names = T, quote=F)
}

colnames(percentage_of_unique_genes) = c("condition","total_number_of_unique_DEGs_in_this_condition", "total_number_of_DEGs_in_this_condition", "percentage_of_unique_genes")
write.table(percentage_of_unique_genes, file = "../specific.comparisons/genes_uniquely_regulated_by_each_condition/table_of_percentage_of_unique_genes.txt", col.names = T, row.names = F, quote=F)


#3. Plot this percentage in a bar chart:
#library(ggplot2); library(grid); library(gridExtra)
color_list = c("lightpink", "darkseagreen2","lightseagreen","lightskyblue", "dodgerblue2", "coral1", "orangered3", "blue", "navy", "purple3")

# plot_list = list()
# for (i in 1:length(condition_name_list)){
#     bacteria = c(condition_name_list[i],condition_name_list[i])
#     names = c("Unique", "All DEGs")
#     numbers = c(percentage_of_unique_genes[i,2], percentage_of_unique_genes[i,3])
#     df = data.frame(bacteria,names,numbers)
#     
#     plot = ggplot(df, aes(x = bacteria, width=10, height=7)) +geom_bar(aes(weight=numbers, fill = rev(names), position = 'stack')) + 
#         scale_y_continuous(paste("condition: ",condition_name_list[i],sep="")) + scale_x_discrete(NULL) + scale_fill_manual(values = c(color_list[i],"ivory")) +
#         theme(axis.title = element_text(face = "bold", size = 30), legend.title = element_text(size = 30), legend.text = element_text(size=30), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
#         theme_bw(base_size = 35)
#     
#     plot_list[[i]] = plot
# }
# 
# grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 5)))
# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# 
# for (n in 1:10){
#     if (n < 6){
#         print(plot_list[[n]], vp = vplayout(1,n))
#     }
#     else {
#         print(plot_list[[n]], vp = vplayout(2,n-5))
#     }
# }

#or
percentage_table = c(percentage_of_unique_genes[,4]);
percentage_list = as.numeric(percentage_table)
percentage_list = percentage_list*100

barplot(percentage_list, main="Percentage of unique DEGs per infection condition", xlab="Infection conditions", ylab = "Percentage of unique DEGs (%)", names.arg=condition_name_list, col= color_list, ylim=c(0,50))






#============================================
#old code from 2/5/2016
#============================================


#1. From the data.table, find the genes that were significant in Staph
on.when.infected.with.staph.12hr = subset(data.table, grepl("Y",data.table[,46]) ); print(dim(on.when.infected.with.staph.12hr)[1])#45,46 (Staph 12hr)
on.only.when.infected.with.staph.12hr = on.when.infected.with.staph.12hr; print(dim(on.only.when.infected.with.staph.12hr)[1])

#2(A). Now with the list of genes that were significant in Staph, only sort genes that have no significant expression in any other time points and any other infection conditions.
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/staph.only")

conditions.other.than.staph.12hr = c(seq(4,44,2),seq(48,64,2))
for (n in conditions.other.than.staph.12hr){
    on.only.when.infected.with.staph.12hr = subset(on.only.when.infected.with.staph.12hr, grepl("N", on.only.when.infected.with.staph.12hr[,n]))
    print(dim(on.only.when.infected.with.staph.12hr)[1])
}
on.only.when.infected.with.staph.12hr = on.only.when.infected.with.staph.12hr[,c(1,2,45)]
on.only.when.infected.with.staph.12hr.sorted = on.only.when.infected.with.staph.12hr[order(on.only.when.infected.with.staph.12hr[,3], decreasing=T),]
write.table(on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] > 0),], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_upregulated.txt", 
            quote=F, col.names= T, row.names=F)
write.table(on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] < 0),], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_downregulated.txt", 
            quote=F, col.names= T, row.names=F)
updegs = on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] > 0),]
write.table(updegs[,1], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_upregulated_gene_name_only.txt", 
            quote=F, col.names= T, row.names=F)

#2/2/2016 -- compare 12hr Staph to 12hr other conditions:
#2(B). Now with the list of genes that were significant in Staph, only pick genes that have no significant expression in any other infection conditions at 12hr.
conditions.other.than.staph.12hr = c(4,10,16,22,28,34,40,48,50,52,54,60) #includes heatkilled bacteria
for (n in conditions.other.than.staph.12hr){
    on.only.when.infected.with.staph.12hr = subset(on.only.when.infected.with.staph.12hr, grepl("N", on.only.when.infected.with.staph.12hr[,n]))
    print(dim(on.only.when.infected.with.staph.12hr)[1])
}
on.only.when.infected.with.staph.12hr = on.only.when.infected.with.staph.12hr[,c(1,2,45)]
on.only.when.infected.with.staph.12hr.sorted = on.only.when.infected.with.staph.12hr[order(on.only.when.infected.with.staph.12hr[,3], decreasing=T),]
head(on.only.when.infected.with.staph.12hr.sorted)

write.table(on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] > 0),], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_at_12hr_upregulated.txt", 
            quote=F, col.names= T, row.names=F)
write.table(on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] < 0),], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_at_12hr_downregulated.txt", 
            quote=F, col.names= T, row.names=F)
updegs = on.only.when.infected.with.staph.12hr.sorted[which(on.only.when.infected.with.staph.12hr.sorted[,3] > 0),]
write.table(updegs[,1], file="list_of_genes_only_activated_by_S.aureus_but_not_by_others_at_12hr_upregulated_gene_name_only.txt", 
            quote=F, col.names= T, row.names=F)


#======What if I do the same with M.luteus?

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data.table <- read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt",header=T)

#1. From the data.table, find the genes that were significant in M.luteus
on.when.infected.with.m.luteus.12hr = subset(data.table, grepl("Y",data.table[,10]) ); print(dim(on.when.infected.with.m.luteus.12hr)[1]) #9,10 (Staph 12hr)
on.only.when.infected.with.m.luteus.12hr = on.when.infected.with.m.luteus.12hr; print(dim(on.only.when.infected.with.m.luteus.12hr)[1])
conditions.other.than.m.luteus.12hr = c(4,16,22,28,34,40,46,48,50,52,54,60) #includes heatkilled bacteria
for (n in conditions.other.than.m.luteus.12hr){
    on.only.when.infected.with.m.luteus.12hr = subset(on.only.when.infected.with.m.luteus.12hr, grepl("N", on.only.when.infected.with.m.luteus.12hr[,n]))
    print(dim(on.only.when.infected.with.m.luteus.12hr)[1])
}
on.only.when.infected.with.m.luteus.12hr = on.only.when.infected.with.m.luteus.12hr[,c(1,2,45)]
on.only.when.infected.with.m.luteus.12hr.sorted = on.only.when.infected.with.m.luteus.12hr[order(on.only.when.infected.with.m.luteus.12hr[,3], decreasing=T),]
head(on.only.when.infected.with.m.luteus.12hr.sorted)

write.table(on.only.when.infected.with.m.luteus.12hr.sorted[which(on.only.when.infected.with.m.luteus.12hr.sorted[,3] > 0),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/staph.only/list_of_genes_only_activated_by_M.luteus_but_not_by_others_at_12hr_upregulated.txt", 
            quote=F, col.names= T, row.names=F)
write.table(on.only.when.infected.with.m.luteus.12hr.sorted[which(on.only.when.infected.with.m.luteus.12hr.sorted[,3] < 0),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/staph.only/list_of_genes_only_activated_by_M.luteus_but_not_by_others_at_12hr_downregulated.txt", 
            quote=F, col.names= T, row.names=F)
updegs = on.only.when.infected.with.m.luteus.12hr.sorted[which(on.only.when.infected.with.m.luteus.12hr.sorted[,3] > 0),]
write.table(updegs[,1], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/staph.only/list_of_genes_only_activated_by_M.luteus_but_not_by_others_at_12hr_upregulated_gene_name_only.txt", 
            quote=F, col.names= T, row.names=F)


