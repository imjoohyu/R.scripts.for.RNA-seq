#Find the recovery genes based on edgeR expression path assignment
#October 25, 2016
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

#Read in data that had already removed genes with EE-EE-EE path.
M.luteus = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_M.luteus_EEs_removed.txt", header=T)
E.coli = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_E.coli_EEs_removed.txt", header=T)
S.mar.type = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_S.mar.type_EEs_removed.txt", header=T)
E.fae.live = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_E.fae.live_EEs_removed.txt", header=T)
P.rett.live = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_P.rett.live_EEs_removed.txt", header=T)
Ecc15 = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_Ecc15_EEs_removed.txt", header=T)
E.fae.heatkilled = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_E.fae.heatkilled_EEs_removed.txt", header=T)
P.rett.heatkilled = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_P.rett.heatkilled_EEs_removed.txt", header=T)
list.of.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, E.fae.heatkilled, P.rett.heatkilled)
list.of.conditions.names <- c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")


#Q1. What are the recovery genes in our dataset?

#Select genes that have EE at the 0h-5.5d comparison (X-X-EE). Do it separately by infection conditions.
possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
recovery_genes = c(); recovery_genes_total = c()
for (i in 1:length(list.of.conditions)){
    data.to.check = as.matrix(as.data.frame(list.of.conditions[i]))
    recovery_genes = c()
    for (j in 1:dim(data.to.check)[1]){
        pattern = as.character(data.to.check[j,3])
        if (pattern %in% possible_patterns){ #if the pattern is one of the 
            recovery_genes = rbind(recovery_genes, data.to.check[j,])
        }
    }
    colnames(recovery_genes) = c("gene_id", "gene_name", paste(list.of.conditions.names[i],".path",sep=""))
    write.table(recovery_genes, file=paste("recovery_genes/recovery_genes_in_",list.of.conditions.names[i],".txt",sep=""), row.names = F, col.names = T, quote=F)
    recovery_genes_total = rbind(recovery_genes_total, recovery_genes)
}

#list of genes that were recovery genes in at least one condition (live and heatkilled)
recovery_genes_total = recovery_genes_total[,c(1,2)]
recovery_genes_total_uniq = unique(recovery_genes_total)
write.table(recovery_genes_total_uniq, file=paste("recovery_genes/recovery_genes_in_at_least_one_condition",".txt",sep=""), row.names = F, col.names = T, quote=F)

#Q2. Do they overlap with the recovery genes from the Schneider paper?
Schneider = read.table("recovery_genes/Schneider_recovery_genes_Flybase_converted.txt", header=T)
recovery_genes = data.frame(recovery_genes_total_uniq)

#overlap by FBgn number
overlap.between.mine.and.Schneider=c()
for (i in 1:dim(recovery_genes)[1]){
    pattern = as.character(recovery_genes[i,1])
    if (pattern %in% Schneider[,2]){
        overlap.between.mine.and.Schneider = rbind(overlap.between.mine.and.Schneider, recovery_genes[i,])
    }
}
dim(overlap.between.mine.and.Schneider) #163
write.table(overlap.between.mine.and.Schneider, file="recovery_genes/overlap_between_flysick_and_Schneider_recovery_genes.txt", row.names = F, col.names = T, quote=F)


#Q3. Do cleanprick DEGs recover? -- find overlap between cleanprick and recovery_genes
#Assign the expression paths for clean prick data
total_de_list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #**Using the data from 2015

#Assign directions
length.of.table <- as.numeric(dim(total_de_list)[1]); width.of.table <- as.numeric(dim(total_de_list)[2]); indicator <- NULL

for (i in 1:length.of.table){ #1, 2, 3, ... 2589
    #cat("The value for k is: ", k, "\n" )
    for (s in seq(3, width.of.table, 2)){ #3, 5, 7, ... 63
        
        if (total_de_list[i,s+1] == "Y"){ #if significant
            indicator <- isTRUE(total_de_list[i,s] > 0) #indicator shows that the direction is positive
            #cat("indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the case is Up-DEG
                total_de_list[i,s] = "Up"
            }
            else { #If the caseis Down-DEG
                total_de_list[i,s] = "Down"
            }
        }
        else { #if not significant
            total_de_list[i,s] = "EE"
        }
        
    }
}

cleanprick = total_de_list[,c(3:8)]

#Put the expression path together and only pick the ones that recover.
cleanprick.expression.path.table = matrix(NA, nrow = dim(cleanprick)[1], ncol= 1)
for (b in 1:dim(cleanprick)[1]){ #for each gene
    subset = cleanprick[b,]
    path = paste(subset[1,1],"-",subset[1,3],"-",subset[1,5], sep="")
    cleanprick.expression.path.table[b,1] = path
}
cleanprick.expression.path.table = cbind(total_de_list[,c(1:2)],cleanprick.expression.path.table)
colnames(cleanprick.expression.path.table) = c("gene_id","gene_name","cleanprick.path")

possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
recovery_genes = c()
for (j in 1:dim(cleanprick.expression.path.table)[1]){
    pattern = as.character(cleanprick.expression.path.table[j,3])
    if (pattern %in% possible_patterns){ #if the pattern is one of the 
        recovery_genes = rbind(recovery_genes, cleanprick.expression.path.table[j,])
    }
}
colnames(recovery_genes) = c("gene_id", "gene_name", "cleanprick.path") #94 genes
write.table(recovery_genes, file=paste("recovery_genes/recovery_genes_in_cleanprick.txt",sep=""), row.names = F, col.names = T, quote=F)


clean.prick.DEGs = data.frame(read.table("recovery_genes/list.of.all.degs.for.cleanprick.txt")) #UC-cleanprick, made in overlap_between_infection_cleanprick_heatkilled.R

#overlap by FBgn number
a=c()
for (i in 1:dim(recovery_genes)[1]){
    pattern = as.character(recovery_genes[i,1])
    if (pattern %in% clean.prick.DEGs[,1]){
        a = rbind(a, recovery_genes[i,])
    }
    else{
        print(pattern)
    }
}
dim(a)
#
