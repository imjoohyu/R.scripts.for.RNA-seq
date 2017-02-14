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


#The other way to do this:
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
        else{
            AMP_data_pattern_checked[j,k] = as.character("No")
        }
    }
}

AMP_data_converted = cbind(AMP_data[,1:2], AMP_data_pattern_checked[,3:11])
colnames(AMP_data_converted) = c("gene_id", "gene_name", "Sterile Wound", "M.lutues", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "E.faecalis heatkilled", "P.rettgeri heatkilled")
write.table(AMP_data_converted, file="recovery_genes/Are_AMPs_recovered.txt", row.names = F, col.names = T, quote =F)


library(gplots)
AMP_data_matrix = as.matrix(AMP_data_converted)
heatmap.2(x=AMP_data_matrix, Rowv= FALSE, Colv = FALSE, dendrogram = "none", cellnote = AMP_data_matrix)




