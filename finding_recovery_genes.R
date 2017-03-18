#Find the recovery genes based on edgeR expression path assignment
#December 10th, 2016
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")


#1. Get the genes with excluding those that were EE-EE-EE in every condition.
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

#Put the expression path together
put_paths_together = function(data){
    full.expression.path.table = matrix(NA, nrow = dim(data)[1], ncol= 9)
    for (b in 1:dim(data)[1]){ #for each gene
        count=1
        for (c in c(3,9,15,21,27,33,39,53,59)) { #for each condition of cleanprick, 6 live conditions and 2 heatkilled conditions 
            start = c; end = c+5
            subset = data[b,start:end]; subset
            path = paste(subset[1,1],"-",subset[1,3],"-",subset[1,5], sep="")
            full.expression.path.table[b, count] = path
            count=count+1
        }
    }
    return(full.expression.path.table) 
}
path_data = put_paths_together(direction_data)
path_data = cbind(direction_data[,c(1:2)], path_data)

#Remove any genes that havent' changed significantly in any of the conditions.
remove_NS = function(data){
    full.expression.path.table.EEs.removed = c()
    for (k in 1:dim(data)[1]){ #1, 2, 3, ... 11911
        for (m in seq(3, 11, 1)){ #1, 4, ..., 9
            indicator <- isTRUE(data[k,m] != "EE-EE-EE") #When the gene is up or down
            if (indicator == T){
                full.expression.path.table.EEs.removed = rbind(full.expression.path.table.EEs.removed, data[k,])
                break
            }
        }
    }
    return(full.expression.path.table.EEs.removed)
}
path_data_NS_removed = remove_NS(path_data) 
colnames(path_data_NS_removed) = c("gene.id","gene.name","clean.prick", "M.luteus", "E.coli", "S.marcescens", "E.faecalis","P.rettgeri","Ecc15", "E.fae.heatkilled","P.rett.heatkilled")


#2. For each condition, count the genes with a path other than EE-EE-EE and count the recovery genes
#Each condition would have different numbers for these two categories.

#Count the genes with a path other than EE-EE-EE
list.of.conditions.names <- c("SterileWound", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")
get_active_DEGs = function(data){
    active_DEGs=matrix(NA, ncol=2, nrow=length(list.of.conditions.names))
    for (i in 1:length(list.of.conditions.names)){ #in a given condition
        data_condition = data[,c(1,2,i+2)]
        #get the total number of active DEGs by removing the genes with the EE-EE-EE path.
        data_condition_active = data_condition[which(data_condition[,3] != "EE-EE-EE"),]
        active_DEGs[i,1] = as.character(list.of.conditions.names[i])
        active_DEGs[i,2] = dim(data_condition_active)[1]
    }
    return(active_DEGs)
}
active_DEGs_table = get_active_DEGs(path_data_NS_removed)

#Select genes that have EE at the 0h-5.5d comparison (X-X-EE). Do it separately by infection conditions.
possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
get_recovery_DEGs = function(data){
    recovery_genes_total = c()
    recovery_DEGs=matrix(NA, ncol=2, nrow=length(list.of.conditions.names))
    
    for (i in 1:length(list.of.conditions.names)){ #in a given condition
        cat("condition: ", list.of.conditions.names[i])
        recovery_genes = c()
        data_condition = data[,c(1,2,i+2)]

        for (j in 1:dim(data_condition)[1]){
            pattern = as.character(data_condition[j,3])
            if (pattern %in% possible_patterns){ #if the pattern is one of the 
                recovery_genes = rbind(recovery_genes, data_condition[j,])
            }
        }
        recovery_DEGs[i,1] = as.character(list.of.conditions.names[i])
        recovery_DEGs[i,2] = as.numeric(dim(recovery_genes)[1])
        
        colnames(recovery_genes) = c("gene_id", "gene_name", paste(list.of.conditions.names[i],".path",sep=""))
        #write.table(recovery_genes, file=paste("recovery_genes/recovery_genes_in_",list.of.conditions.names[i],".txt",sep=""), row.names = F, col.names = T, quote=F)
        
        #Put the recovery genes across the conditions together
        colnames(recovery_genes)=c("gene_id","gene_name","gene_path")
        recovery_genes_total = rbind(recovery_genes_total, recovery_genes)
    }
    
    #list of genes that were recovery genes in at least one condition (live and heatkilled)
    recovery_genes_total = recovery_genes_total[,c(1,2)]
    recovery_genes_total_uniq = unique(recovery_genes_total)
    #write.table(recovery_genes_total_uniq, file="recovery_genes/recovery_genes_in_at_least_one_condition.txt", row.names = F, col.names = T, quote=F)

    colnames(recovery_DEGs) = c("condition","number_of_recovery_genes")
    return(recovery_DEGs)
}
recovery_DEGs_table = get_recovery_DEGs(path_data_NS_removed)




#3. Research questions:
#3A. What percentage of the active DEGs recovery genes?
perc_table = data.frame(cbind(as.numeric(as.character(recovery_DEGs_table[,2])), as.numeric(as.character(active_DEGs_table[,2]))))
perc_table[,3] = round((perc_table[,1]/perc_table[,2])*100,2)
perc_table = cbind(active_DEGs_table[,1],perc_table)
colnames(perc_table) = c("condition","number_of_recovery_genes","number_of_DEGs","percentage_of_recovery_genes")

#Add the entry on genes that recovered in at least one condition & percentage
total_recovery_genes = read.table("recovery_genes/recovery_genes_in_at_least_one_condition.txt", header=T)
perc_table$condition = as.character(perc_table$condition)
perc_table[10,1] = as.character("At_least_one"); perc_table[10,2]= dim(total_recovery_genes)[1]
perc_table[10,3] = dim(path_data_NS_removed)[1]; perc_table[10,4] = round((perc_table[10,2]/perc_table[10,3])*100,2)

#Create a stacked percentage bar plot
library(ggplot2); library(scales); library(reshape2)
perc_table[,5] = perc_table[,3] - perc_table[,2]
colnames(perc_table)[5] = c("number_of_DEGs_not_recovered")

datm <- melt(cbind(perc_table[,c(2,5)], condition = perc_table$condition), id.vars = c('condition'))
pct = as.character(perc_table$percentage_of_recovery_genes)
positions = c("SterileWound", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "E.fae.heatkilled", "P.rett.live", "P.rett.heatkilled", "Ecc15","At_least_one")
ggplot(datm,aes(x = condition,y = value, fill = variable)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + guides(fill=FALSE) + scale_x_discrete(limits = positions) + scale_fill_manual(values = c("skyblue", "grey"))

#The other way to plot this graph with different colors representing different bacterial infections:
positions = c("SterileWound", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "E.fae.heatkilled", "P.rett.live", "P.rett.heatkilled", "Ecc15","At_least_one")
ggplot(data=perc_table, aes(x=condition, y=percentage_of_recovery_genes, fill=condition)) + geom_bar(stat="identity") + scale_fill_manual(values = c("lightgrey", "darkseagreen2", "dark salmon", "coral1", "lightskyblue", "lightpink", "steelblue1", "dodgerblue2", "lightseagreen","azure4"))+ scale_x_discrete(limits = positions) +  scale_y_continuous(limits=c(0,100)) + labs(y = "Percentage (%) of recovery genes") + theme_bw(base_size=14) + guides(fill=FALSE) + geom_text(size=6, aes(label=paste(percentage_of_recovery_genes,"%")))


#+ geom_text(aes(x=paste0(pct,"%")), size=4) -- add % to the graph



#No need to check this separately (12/11/2016)
# #3B. What percentage of the clean prick genes recover?
# clean.prick.DEGs = data.frame(read.table("recovery_genes/list.of.all.degs.for.cleanprick.txt")) #UC-cleanprick, made in overlap_between_infection_cleanprick_heatkilled.R
# clean.prick.recovery.genes = read.table("recovery_genes/recovery_genes_in_SterileWound.txt",header=T)
# 
# #overlap by FBgn number
# a=c()
# for (i in 1:dim(clean.prick.recovery.genes)[1]){
#     pattern = as.character(clean.prick.recovery.genes[i,1])
#     if (pattern %in% clean.prick.DEGs[,1]){
#         a = rbind(a, clean.prick.recovery.genes[i,])
#     }
#     else{
#         print(pattern)
#     }
# }
# dim(a) #everything overlaps


#3B. What percentage of recovery genes overlap with the Schneider recovery genes?
Schneider = read.table("recovery_genes/Schneider_recovery_genes_Flybase_converted.txt", header=T)

#overlap by FBgn number
overlap.between.mine.and.Schneider=c()
for (i in 1:dim(total_recovery_genes)[1]){
    pattern = as.character(total_recovery_genes[i,1])
    if (pattern %in% Schneider[,2]){
        overlap.between.mine.and.Schneider = rbind(overlap.between.mine.and.Schneider, total_recovery_genes[i,])
    }
}
dim(overlap.between.mine.and.Schneider) #167
write.table(overlap.between.mine.and.Schneider, file="recovery_genes/overlap_between_flysick_and_Schneider_recovery_genes.txt", row.names = F, col.names = T, quote=F)



#3C. Are recovery genes specific to a condition or are ‘common’ and thus occur in many conditions?
get_freq_of_recovery_genes = function(data, possible_patterns){
    specificity_table = matrix(NA, nrow=dim(data)[1], ncol =3); count=0; count_total=0
    
    for (i in 1:dim(data)[1]){
        gene = data[i,c(3:dim(data)[2])]; gene.t = t(gene)
        
        for (j in 1:length(possible_patterns)){
            pattern = possible_patterns[j]
            for (k in 1:dim(data)[2]){
                if (pattern %in% gene[1,k]){#if there is at least one condition with the said pattern
                    count = count + 1
                }
            }
            count_total = count_total + count
            count=0
        }
        
        specificity_table[i,1] = as.character(data[i,1])
        specificity_table[i,2] = as.character(data[i,2])
        specificity_table[i,3] = count_total
        count=0; count_total=0
    }
    return(specificity_table)
}

#For all 9 conditions
possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
specificity_table = get_freq_of_recovery_genes(path_data_NS_removed, possible_patterns)

#For live infections only
path_data_NS_removed_live_only = path_data_NS_removed[,c(1,2,4,5,6,7,8,9)]
possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
specificity_table_live_infection = get_freq_of_recovery_genes(path_data_NS_removed_live_only, possible_patterns)

#Only get the genes that had a recovery pattern in at least one condition
specificity_table_removing_zero = specificity_table[which(specificity_table[,3] != 0),]
specificity_table_live_infection_removing_zero = specificity_table_live_infection[which(specificity_table_live_infection[,3] != 0),]

#Plot the histogram
barplot(table(specificity_table_removing_zero[,3]), ylim=c(0,1000), main= "All conditions", xlab = "Number of conditions", ylab = "Number of recovery genes")
barplot(table(specificity_table_live_infection_removing_zero[,3]), ylim=c(0,1000), main= "Live conditions only", xlab = "Number of conditions", ylab = "Number of recovery genes")


#What if splitting up-then-recovered genes (Up-EE-EE, EE-Up-EE, Up-Up-EE) and down-then-recovered genes (Down-EE-EE, EE-Down-EE, Down-Down-EE)?
possible_patterns = c("Up-Up-EE","Up-EE-EE","EE-Up-EE")
up_recovery_genes = get_freq_of_recovery_genes(path_data_NS_removed, possible_patterns)
up_recovery_genes_removing_zero = up_recovery_genes[which(up_recovery_genes[,3] != 0),]
possible_patterns = c("EE-Down-EE","Down-EE-EE", "Down-Down-EE")
down_recovery_genes = get_freq_of_recovery_genes(path_data_NS_removed, possible_patterns)
down_recovery_genes_removing_zero = down_recovery_genes[which(down_recovery_genes[,3] != 0),]

barplot(table(up_recovery_genes_removing_zero[,3]), ylim=c(0,500), main= "Genes upregulated then recovered", xlab = "Number of conditions", ylab = "Number of recovery genes")
barplot(table(down_recovery_genes_removing_zero[,3]), ylim=c(0,500), main= "Genes downregulated then recovered", xlab = "Number of conditions", ylab = "Number of recovery genes")


#Do some genes that showed an up-then-recovered pattern in some conditions show the opposite pattern in other conditions? Yes
path_with_freq = cbind(path_data_NS_removed,specificity_table[,3],specificity_table_live_infection[,3],up_recovery_genes[,3], down_recovery_genes[,3])
colnames(path_with_freq)[12:15] = c("num_of_rec_cond_from_all_9","num_of_rec_cond_from_live_6","up_recovery_genes_from_all_9","down_recovery_gene_from_all_9")
path_with_freq_removing_zero = path_with_freq[which(path_with_freq$num_of_rec_cond_from_all_9 != 0),]

#A few ways to identify these genes:
genes_with_dual_recovery_pattern = path_with_freq[which(path_with_freq$up_recovery_genes_from_all_9 ==1 & path_with_freq$down_recovery_gene_from_all_9 == 1),]
genes_with_dual_recovery_pattern = path_with_freq[which(path_with_freq$up_recovery_genes_from_all_9 ==2 & path_with_freq$down_recovery_gene_from_all_9 == 1),]
genes_with_dual_recovery_pattern = path_with_freq[which(path_with_freq$up_recovery_genes_from_all_9 ==1 & path_with_freq$down_recovery_gene_from_all_9 == 2),]
genes_with_dual_recovery_pattern = path_with_freq[which(path_with_freq$up_recovery_genes_from_all_9 ==2 & path_with_freq$down_recovery_gene_from_all_9 == 2),]


#What if I remove the uniquely expressed genes from the recovery genes bin 1? (1/13/2017)
unique_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/genes_uniquely_regulated_by_each_condition/unique_genes_ignoring_time.txt", header=T)
unique_genes_simplified = unique_genes[,c(1:2)]
unique_genes_simplified_uniquified = unique(unique_genes_simplified)

#Remove uniquely expressed genes from the recovery genes that only recovered in one condition (bin 1)
remove_express_genes_from_bin1 = function(data){
    bin1 = data[which(data[,3] == "1"),]
    the_rest = data[which(data[,3] != "1"),]
    for (i in 1:dim(unique_genes_simplified_uniquified)[1]){
        gene_to_remove = unique_genes_simplified_uniquified[i,1]
        bin1 = bin1[which(bin1[,1] != gene_to_remove),]
    }
    data_with_new_bin1 = rbind(bin1, the_rest)
    return(barplot(table(data_with_new_bin1[,3]), ylim=c(0,500), xlab = "Number of conditions", ylab = "Number of recovery genes"))
    
}
par(mfrow = c(2,2))
specificity_table_all_conditions_unique_genes_removed = remove_express_genes_from_bin1(specificity_table_removing_zero) #all conditions
specificity_table_all_conditions_unique_genes_removed = remove_express_genes_from_bin1(up_recovery_genes_removing_zero) #genes upregulated than recovered
specificity_table_live_infection_unique_genes_removed = remove_express_genes_from_bin1(specificity_table_live_infection_removing_zero) #live conditions only
specificity_table_all_conditions_unique_genes_removed = remove_express_genes_from_bin1(down_recovery_genes_removing_zero) #genes downregulated then recovered



#3D. What are the genes that were Up-Down-EE or Down-Up-EE?
possible_patterns = c("Up-Down-EE","Down-Up-EE")
odd_recovery_genes = get_freq_of_recovery_genes(path_data_NS_removed, possible_patterns)
table(odd_recovery_genes[,3]) #0
odd_recovery_genes = get_freq_of_recovery_genes(path_data, possible_patterns)
table(odd_recovery_genes[,3]) #0


#3E. What’s the difference between the GO term of recovery vs non-recovery genes?
#i) Genes that recovered in at least one condition vs genes that were active in at least one condition but were not recovered
total_recovery_genes = read.table("recovery_genes/recovery_genes_in_at_least_one_condition.txt", header=T)
data_to_trim = path_data_NS_removed
non_recovery_genes=c()
for (m in 1:dim(total_recovery_genes)[1]){
    gene_to_look_at = as.character(total_recovery_genes[m,1])
    data_to_trim = data_to_trim[!(data_to_trim$gene.id == gene_to_look_at),]
    non_recovery_genes = data_to_trim
}
write.table(non_recovery_genes[,1:2], file="recovery_genes/non_recovery_genes.txt", row.names = F, col.names = T, quote=F)

#ii) Genes that recover in multiple conditions (the right side of the hisgram/barplot) = so-called “super recovery genes"
super_core_genes_7_plus = path_with_freq_removing_zero[which(path_with_freq_removing_zero$num_of_rec_cond_from_all_9 == 9 | path_with_freq_removing_zero$num_of_rec_cond_from_all_9 == 8 | path_with_freq_removing_zero$num_of_rec_cond_from_all_9 == 7),] #35 genes
super_core_genes_5_plus_live_only = path_with_freq_removing_zero[which(path_with_freq_removing_zero$num_of_rec_cond_from_live_6 == 6 | path_with_freq_removing_zero$num_of_rec_cond_from_live_6 == 5),] #79 genes
write.table(super_core_genes_7_plus, file="recovery_genes/genes_that_recovered_in_7_or_more_all_conditions.txt", row.names = F, col.names = T, quote=F)
write.table(super_core_genes_5_plus_live_only , file="recovery_genes/genes_that_recovered_in_5_or_more_live_conditions.txt", row.names = F, col.names = T, quote=F)


#3F. What happens to the AMPs? (2/13/2017)
input_data= path_data_NS_removed
AMPs = c("Dro", "Drs", "CecA1", "CecA2", "CecB", "CecC", "Dpt", "DptB", "AttA", "AttB", "Def", "Mtk")
possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
no_change = c("EE-EE-EE")
AMP_data = c()
for (i in 1:length(AMPs)){
    AMP_name = as.character(AMPs[i])
    AMP_data = rbind(AMP_data, input_data[which(input_data$gene.name == AMP_name),])
}
AMP_data_pattern_checked = matrix(NA, nrow=dim(AMP_data)[1], ncol=dim(AMP_data)[2])
for (j in 1:length(AMPs)){
    for (k in 3:11){
        if (as.character(AMP_data[j,k]) %in% possible_patterns){ #if recovered/returned
            AMP_data_pattern_checked[j,k] = as.character("Yes")
        }
        else if (as.character(genes_to_test_data[j,k]) == no_change) { #if no change
            genes_to_test_data_checked[j,k] = as.character("No_change")
        }
        else{
            AMP_data_pattern_checked[j,k] = as.character("No")
        }
    }
}

AMP_data_converted = cbind(AMP_data[,1:2], AMP_data_pattern_checked[,3:11])
colnames(AMP_data_converted) = c("gene_id", "gene_name", "Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "E.faecalis heatkilled", "P.rettgeri heatkilled")
write.table(AMP_data_converted, file="recovery_genes/Are_AMPs_recovered.txt", row.names = F, col.names = T, quote =F)


#What about other functions? (2/27/2017, Updated on 3/18/2017)
#generic function (input: AmiGo2 output):
library(ggplot2); library(reshape2); library(klaR); library(RColorBrewer)
check_recovery_status_by_category = function(address_to_test_category, include) {
    
    if(include == "Yes"){
        input_data = path_data #includes genes that did not change across all conditions
        colnames(input_data) = colnames(path_data_NS_removed)
    }
    else if(include == "No"){
        input_data= path_data_NS_removed #excludes genes that did not change across all conditions
    }
    
    genes_to_test = read.table(address_to_test_category, header=F, sep="\t")
    genes_to_test[,1] = sapply(genes_to_test[,1], function(x) {sub("FB:","",x)} )
    genes_to_test = unique(genes_to_test)
    
    possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
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
            if (as.character(genes_to_test_data[j,k]) %in% possible_patterns){ #if recovered/returned
                genes_to_test_data_checked[j,k] = as.character("Yes")
            }
            else if (as.character(genes_to_test_data[j,k]) == no_change) { #if no change
                genes_to_test_data_checked[j,k] = as.character("No_change")
            }
            else{
                genes_to_test_data_checked[j,k] = as.character("No")
            }
        }
    }
    data_converted = cbind(genes_to_test_data[,1:2], genes_to_test_data_checked[,3:dim(genes_to_test_data_checked)[2]])
    colnames(data_converted) = c("gene_id", "gene_name", "Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "E.faecalis heatkilled", "P.rettgeri heatkilled")
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
plot_the_recovery_status = function(recovery_status_table, sort_by_gene, clustering_info, condition_name){
    data = recovery_status_table
    
    if (sort_by_gene == "Yes"){
        condition_order = c("Sterile Wound","M.luteus","E.coli","S.marcescens Type","E.faecalis live","P.rettgeri live", "Ecc15", "E.faecalis heatkilled","P.rettgeri heatkilled")
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
    else{
        condition_order = c("Sterile Wound","M.luteus","E.coli","S.marcescens Type","E.faecalis live","P.rettgeri live", "Ecc15", "E.faecalis heatkilled","P.rettgeri heatkilled")
        condition_order_with_clustering_info = cbind(condition_order, clustering_info)
        condition_order_with_clustering_info = condition_order_with_clustering_info[order(condition_order_with_clustering_info[,2], decreasing = T),]
        condition_order = condition_order_with_clustering_info[,1]
        
        data$gene_id = droplevels(data$gene_id); data$gene_name = droplevels(data$gene_name)
        data_rearranged <- melt(data, id = c('gene_id', 'gene_name'))
        data_rearranged$variable = factor(data_rearranged$variable, levels= condition_order)
        colnames(data_rearranged) = c("gene_id", "gene_name", "condition", "recovery_pattern")
    }
    
    plot = ggplot(data_rearranged, aes(condition, gene_name, recovery_pattern)) + geom_tile(aes(fill = recovery_pattern), colour = "white") + scale_fill_manual(values=c("#d8b365", "lightgrey", "#5ab4ac")) + ggtitle(paste("Recovery pattern for ", condition_name, sep="") ) + theme(text=element_text(size=15))
    return(plot)
}
visualize_patterns = function(address_to_test_category, include, sort_by_gene, number_of_cluster, condition_name){
    response_data = check_recovery_status_by_category(address_to_test_category, include)
    clustering_info = cluster_the_genes_by_pattern(response_data, sort_by_gene, number_of_cluster)
    plot_the_recovery_status(response_data, sort_by_gene, clustering_info, condition_name)
} #all put together


#1. Response to Wounding -- 38 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(38)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_wounding_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "response to wounding")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_wounding_genes_based_on_AmiGO.txt", "No", "No", 4, "response to wounding")


#2. Cuticle development -- 32 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(32)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/cuticle_development_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "cuticle development")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/cuticle_development_genes_based_on_AmiGO.txt", "No", "No", 4, "cuticle development")


#3. Response to oxidative stress -- 24 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(24)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_oxidative_stress_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "oxidative stress")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_oxidative_stress_genes_based_on_AmiGO.txt", "No", "No", 4, "oxidative stress")


#4. Humoral immune response/AMP production(lower-level category)
#Humoral immune response -- 53 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(53)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/humoral_immune_response_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "humoral immune response")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/humoral_immune_response_genes_based_on_AmiGO.txt", "No", "No", 4, "humoral immune response")

#AMP production -- 16 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(16)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/AMP_production_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "AMP")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/AMP_production_genes_based_on_AmiGO.txt", "No", "No", 4, "AMP")


#5. Metal ion transport and homeostasis -- 34 genes
#Cluster the patterns by genes
cluster_number = ceiling(sqrt(34)/2) #determine the number of clusters for clustering by genes
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/metal_ion_merged_genes_based_on_AmiGO.txt", "No", "Yes", cluster_number, "metal ion transport and homeostasis")
#Cluster the patterns by conditions
visualize_patterns("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/metal_ion_merged_genes_based_on_AmiGO.txt", "No", "No", 4, "metal ion transport and homeostasis")








#---
#6. Imd pathway between E.coli (benign) and Ecc15/P.rettgeri (chronic)
imd_benign_chronic = check_recovery_status_by_category("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/response_to_Gram-negative_potential_Imd_genes_based_on_AmiGO.txt", "Yes") 
condition_order = c("E.coli","Ecc15","P.rettgeri live"); condition_name = "Imd pathway"
imd_benign_chronic = imd_benign_chronic[,c(1,2,5,9,8)]
imd_benign_chronic$gene_id = droplevels(imd_benign_chronic$gene_id); imd_benign_chronic$gene_name = droplevels(imd_benign_chronic$gene_name)
data_rearranged <- melt(imd_benign_chronic, id = c('gene_id', 'gene_name'))
data_rearranged$variable = factor(data_rearranged$variable, levels= condition_order)
colnames(data_rearranged) = c("gene_id", "gene_name", "condition", "recovery_pattern")
plot = ggplot(data_rearranged, aes(condition, gene_name, recovery_pattern)) + geom_tile(aes(fill = recovery_pattern), colour = "white") + scale_fill_manual(values=c("#d8b365", "lightgrey", "#5ab4ac")) + ggtitle(paste("Recovery pattern for ", condition_name, sep="") ) + theme(text=element_text(size=15))
plot

print_variable_name <- function(x) {
    deparse(substitute(x))
}





##
data = recovery_status_table
data$gene_id = droplevels(data$gene_id); data$gene_name = droplevels(data$gene_name)