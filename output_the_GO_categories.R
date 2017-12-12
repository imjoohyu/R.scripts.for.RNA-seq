#Processing GO terms
#February 18th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry

#script
#Output the GO categories associated with each gene
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


