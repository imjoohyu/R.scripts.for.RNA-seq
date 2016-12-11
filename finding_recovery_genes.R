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
        write.table(recovery_genes, file=paste("recovery_genes/recovery_genes_in_",list.of.conditions.names[i],".txt",sep=""), row.names = F, col.names = T, quote=F)
        
        #Put the recovery genes across the conditions together
        colnames(recovery_genes)=c("gene_id","gene_name","gene_path")
        recovery_genes_total = rbind(recovery_genes_total, recovery_genes)
    }
    
    #list of genes that were recovery genes in at least one condition (live and heatkilled)
    recovery_genes_total = recovery_genes_total[,c(1,2)]
    recovery_genes_total_uniq = unique(recovery_genes_total)
    write.table(recovery_genes_total_uniq, file="recovery_genes/recovery_genes_in_at_least_one_condition.txt", row.names = F, col.names = T, quote=F)

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
library(ggplot2); library(scales)
perc_table[,5] = perc_table[,3] - perc_table[,2]
colnames(perc_table)[5] = c("number_of_DEGs_not_recovered")

datm <- melt(cbind(perc_table[,c(2,5)], condition = perc_table$condition), id.vars = c('condition'))
pct = as.character(perc_table$percentage_of_recovery_genes)
positions = c("SterileWound", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "E.fae.heatkilled", "P.rett.live", "P.rett.heatkilled", "Ecc15","At_least_one")
ggplot(datm,aes(x = condition,y = value, fill = variable)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + guides(fill=FALSE) + scale_x_discrete(limits = positions) + scale_fill_manual(values = c("skyblue", "grey"))

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







