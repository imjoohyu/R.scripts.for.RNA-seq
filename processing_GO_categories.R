#Processing GO terms
#February 18th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry


#Tool 1: Output the GO categories associated with each gene
get_GO_terms_for_each_gene = function(GO_file_address, list_of_genes_file_address){
    GO_output = read.table(GO_file_address, header=T, sep="\t")
    list_of_genes = read.table(list_of_genes_file_address, header=T)
    GO_output$Genes = as.character(GO_output$Genes)
    
    table_of_gene_with_GO = c()
    #colnames(table_of_gene_with_GO) = c("gene_id", "gene_name", "num.of.sig.infection.conditions", "GO_categories")
    for (i in 1:dim(list_of_genes)[1]){
        gene_id = as.character(list_of_genes[i,1])
        gene_id = sub("FBgn", "FBGN", gene_id)
        GO_categories_table = GO_output[grep(gene_id, GO_output$Genes),]
        if (dim(GO_categories_table)[1] > 0){#if you can find this gene
            GO_categories = toString(as.character(GO_categories_table$Term))
            table_of_gene_with_GO = rbind(table_of_gene_with_GO, cbind(list_of_genes[i,], GO_categories))
        }
        else{
            GO_categories = "NA"
            table_of_gene_with_GO = rbind(table_of_gene_with_GO, cbind(list_of_genes[i,], GO_categories))
        }
    }
    return(table_of_gene_with_GO)
}

#1. Output the GO categories associated with each core upregulated gene
core_upregulated_genes_with_GO_terms = get_GO_terms_for_each_gene("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_GO.txt", "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt")
write.table(core_upregulated_genes_with_GO_terms, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list_matched_with_GOterms.txt", quote=F, col.names = T, row.names = F, sep="\t")

#Output 
a = get_GO_terms_for_each_gene("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others_simplified_GO_DAVID.txt", "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others.txt")


#Tool 2: Output the GO data sorted by percentage and p-value
simplify_the_GO_terms = function(directory_address, direction){
    setwd(directory_address)
    condition_list = c("clean.prick", "M.luteus", "E.coli", "S.mar.type", "Ecc15", "P.rett.live", "E.fae.live", "E.fae.heatkilled", "P.rett.heatkilled","S.aureus", "P.sneebia","S.mar.Db11","P.ento")
    time_list = c("12hr","36hr","5.5d")
    GO_table = c()
    
    for (i in 1:length(condition_list)){
        if (i < 10){ #conditions with three time points
            for (j in 1:length(time_list)){
                cat("condition: ", condition_list[i], ", time: ", time_list[j], "\n")
                input_file = read.table(paste("GO_",condition_list[i],".",time_list[j],"_",direction,".txt",sep=""), header=T, sep="\t", quote="")
                colnames(input_file) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
                input_file_sorted = input_file[which(input_file$Category == 'GOTERM_BP_FAT' | input_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
                #input_file_sorted = input_file_sorted[which(input_file_sorted$PValue < 0.01),] #get rid of entries with EASE p-value less than 0.01
                input_file_sorted = input_file_sorted[order(input_file_sorted$Percentage,decreasing = T),] #sort by percentage -- most common GO terms
                write.table(input_file_sorted, file=paste("GO_",condition_list[i],".",time_list[j],"_",direction,"_most_common_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
                input_file_sorted = input_file_sorted[order(input_file_sorted$PValue,decreasing = F),] #sort by p-value -- most enriched GO terms
                write.table(input_file_sorted, file=paste("GO_",condition_list[i],".",time_list[j],"_",direction,"_most_enriched_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
                
            }
        }
        else { #conditions with one time point
            cat("condition: ", condition_list[i], ", time: 12hr", "\n")
            input_file = read.table(paste("GO_",condition_list[i],".12hr_",direction,".txt", sep=""), header=T, sep="\t", quote="")
            colnames(input_file) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
            input_file_sorted = input_file[which(input_file$Category == 'GOTERM_BP_FAT' | input_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
            #input_file_sorted = input_file_sorted[which(input_file_sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
            input_file_sorted = input_file_sorted[order(input_file_sorted$Percentage,decreasing = T),] #sort by percentage -- most common GO terms
            write.table(input_file_sorted, file=paste("GO_",condition_list[i],".12hr_",direction,"_most_common_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
            input_file_sorted = input_file_sorted[order(input_file_sorted$PValue,decreasing = F),] #sort by p-value -- most enriched GO terms
            write.table(input_file_sorted, file=paste("GO_",condition_list[i],".12hr_",direction,"_most_enriched_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
        }
    }
} #DEGs

#1. Get the GO table for upregulated genes
simplify_the_GO_terms("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_upregulated_DE", "up")
#2. Get the GO table for downregulated genes
simplify_the_GO_terms("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_downregulated_DE", "down")



#Tool 3: Get the GO terms and the genes
condition_list = c("clean.prick", "M.luteus", "E.coli", "S.mar.type", "Ecc15", "P.rett.live", "E.fae.live", "E.fae.heatkilled", "P.rett.heatkilled","S.aureus", "P.sneebia","S.mar.Db11","P.ento")
time_list = c("12hr","36hr","5.5d")
#direction = "up"; setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_upregulated_DE")
direction = "down"; setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_downregulated_DE")

i = 13; j = 1
cat("condition: ", condition_list[i], ", time: ", time_list[j], "\n")
input_file = read.table(paste("GO_",condition_list[i],".",time_list[j],"_",direction,".txt",sep=""), header=T, sep="\t", quote="")
colnames(input_file) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
input_file_sorted = input_file[which(input_file$Category == 'GOTERM_BP_FAT' | input_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
input_file_sorted = input_file_sorted[order(input_file_sorted$Percentage,decreasing = T),] #sort by percentage -- most common GO terms
input_file_sorted = input_file_sorted[order(input_file_sorted$PValue,decreasing = F),] #sort by p-value -- most enriched GO terms

write.table(input_file_sorted[,c(2,6)], file="downregulated_sample_table.txt", quote=F, col.names = T, row.names = F, sep="\t")





###not done as of 2/26/2017
#Tool 3:Pick the most interesting GO terms in each condition
Get_GO_terms_for_each_condition = function(directory_address, direction){
    setwd(directory_address)
    condition_list = c("clean.prick", "M.luteus", "E.coli", "S.mar.type", "Ecc15", "P.rett.live", "E.fae.live", "E.fae.heatkilled", "P.rett.heatkilled","S.aureus", "P.sneebia","S.mar.Db11","P.ento")
    time_list = c("12hr","36hr","5.5d")
    GO_table = c()
    
    for (i in 1:length(condition_list)){
        if (i < 10){ #conditions with three time points
            for (j in 1:length(time_list)){
                cat("condition: ", condition_list[i], ", time: ", time_list[j], "\n")
                input_file = read.table(paste("GO_",condition_list[i],".",time_list[j],"_",direction,".txt",sep=""), header=T, sep="\t", quote="")
                colnames(input_file) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
                input_file_sorted = input_file[which(input_file$Category == 'GOTERM_BP_FAT' | input_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
                #input_file_sorted = input_file_sorted[which(input_file_sorted$PValue < 0.01),] #get rid of entries with EASE p-value less than 0.01
                input_file_sorted = input_file_sorted[order(input_file_sorted$Percentage,decreasing = T),] #sort by percentage -- most common GO terms
                write.table(input_file_sorted, file=paste("GO_",condition_list[i],".",time_list[j],"_",direction,"_most_common_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
                input_file_sorted = input_file_sorted[order(input_file_sorted$PValue,decreasing = F),] #sort by p-value -- most enriched GO terms
                write.table(input_file_sorted, file=paste("GO_",condition_list[i],".",time_list[j],"_",direction,"_most_enriched_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
                
                #                 if (dim(input_file_sorted)[1] == 0){
                #                     data_to_add = data.frame(list(paste(condition_list[i],"_",time_list[j],sep=""), "NA", "NA", "NA"))
                #                     colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                #                     GO_table = rbind(GO_table, data_to_add)
                #                 } #if no data is available
                #                 else{
                #                     data_to_add = data.frame(list(paste(condition_list[i],"_",time_list[j],sep=""), input_file_sorted[,2], input_file_sorted[,4], input_file_sorted[,6]))
                #                     colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                #                     GO_table = rbind(GO_table, data_to_add)
                #                }
                
                #print(input_file_sorted[1:4,1:5]) #is this sorted?
            }
        }
        else { #conditions with one time point
            cat("condition: ", condition_list[i], ", time: 12hr", "\n")
            input_file = read.table(paste("GO_",condition_list[i],".12hr_",direction,".txt", sep=""), header=T, sep="\t", quote="")
            colnames(input_file) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
            input_file_sorted = input_file[which(input_file$Category == 'GOTERM_BP_FAT' | input_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
            #input_file_sorted = input_file_sorted[which(input_file_sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
            input_file_sorted = input_file_sorted[order(input_file_sorted$Percentage,decreasing = T),] #sort by percentage -- most common GO terms
            write.table(input_file_sorted, file=paste("GO_",condition_list[i],".12hr_",direction,"_most_common_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
            input_file_sorted = input_file_sorted[order(input_file_sorted$PValue,decreasing = F),] #sort by p-value -- most enriched GO terms
            write.table(input_file_sorted, file=paste("GO_",condition_list[i],".12hr_",direction,"_most_enriched_GO_terms.txt",sep=""), quote=F, row.names = F, col.names = T) #saved the significant cases
            
            #             if (dim(input_file_sorted)[1] == 0){ #if no data is available
            #                 data_to_add = data.frame(list(paste(condition_list[i],"_12hr", sep=""), "NA", "NA", "NA"))
            #                 colnames(data_to_add) = c("condition", "term", "percentage", "genes")
            #                 GO_table = rbind(GO_table, data_to_add)
            #             }
            #             else{
            #                 
            #                 data_to_add = data.frame(list(paste(condition_list[i],"_12hr", sep=""), input_file_sorted[,2], input_file_sorted[,4], input_file_sorted[,6]))
            #                 colnames(data_to_add) = c("condition", "term", "percentage", "genes")
            #                 GO_table = rbind(GO_table, data_to_add)
            #             }
            #print(input_file_sorted[1:4,1:5])
        }
    }
    #return(GO_table)
}
#code taking in previously sorted documents:
Get_GO_terms_for_each_condition = function(directory_address, direction){
    setwd(directory_address)
    condition_list = c("clean.prick", "M.luteus", "E.coli", "S.mar.type", "Ecc15", "P.rett.live", "E.fae.live", "E.fae.heatkilled", "P.rett.heatkilled","S.aureus", "P.sneebia","S.mar.Db11","P.ento")
    time_list = c("12hr","36hr","5.5d")
    GO_table = c()
    
    for (i in 1:length(condition_list)){
        if (i < 10){ #conditions with three time points
            for (j in 1:length(time_list)){
                GO_file = read.table(paste("GO_full_",condition_list[i],".",time_list[j],"_",direction,"_sorted.txt",sep=""), header=T, sep="\t")
                KEGG_file = read.table(paste("KEGG_GEA_full_",condition_list[i],".",time_list[j],"_",direction,"_sorted.txt",sep=""), header=T, sep="\t")
                if (dim(KEGG_file)[1] > 0){
                    GO_file = rbind(GO_file, KEGG_file) #add KEGG information if available
                }
                GO_file_sorted = GO_file[which(GO_file$Category == 'GOTERM_BP_FAT' | GO_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
                if (dim(GO_file_sorted)[1] == 0){
                    data_to_add = data.frame(list(paste(condition_list[i],"_",time_list[j],sep=""), "NA", "NA", "NA"))
                    colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                    GO_table = rbind(GO_table, data_to_add)
                }
                else{
                    GO_file_sorted =  GO_file_sorted[which(GO_file_sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
                    data_to_add = data.frame(list(paste(condition_list[i],"_",time_list[j],sep=""), GO_file_sorted[,2], GO_file_sorted[,4], GO_file_sorted[,6]))
                    colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                    GO_table = rbind(GO_table, data_to_add)
                }
            }
        }
        else { #conditions with one time point
            GO_file = read.table(paste("GO_full_",condition_list[i],".12hr_",direction,"_sorted.txt", sep=""), header=T, sep="\t")
            KEGG_file = read.table(paste("KEGG_GEA_full_",condition_list[i],".12hr_",direction,"_sorted.txt",sep=""), header=T, sep="\t")
            if (dim(KEGG_file)[1] > 0){
                GO_file = rbind(GO_file, KEGG_file) #add KEGG information if available
            }
            GO_file_sorted = GO_file[which(GO_file$Category == 'GOTERM_BP_FAT' | GO_file$Category == 'KEGG_PATHWAY'),] #only pick biological function
            if (dim(GO_file_sorted)[1] == 0){
                data_to_add = data.frame(list(paste(condition_list[i],"_12hr", sep=""), "NA", "NA", "NA"))
                colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                GO_table = rbind(GO_table, data_to_add)
            }
            else{
                GO_file_sorted =  GO_file_sorted[which(GO_file_sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
                data_to_add = data.frame(list(paste(condition_list[i],"_12hr", sep=""), GO_file_sorted[,2], GO_file_sorted[,4], GO_file_sorted[,6]))
                colnames(data_to_add) = c("condition", "term", "percentage", "genes")
                GO_table = rbind(GO_table, data_to_add)
            }
        }
    }
    return(GO_table)
}
