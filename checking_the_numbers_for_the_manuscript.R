#Checking the numbers for the manuscript
#Febraury 12, 2017
#Joo Hyun Im (ji72)


rm(list=ls(all=TRUE)) #delete any previous entry

#Intriguingly, we also observed a large difference in the number of differentially expressed genes between benign bacteria.


count_DEGs_in_each_condition_without_overlap = function(bacteria_name){
    bacteria1 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_upregulated_DE_genes_in_", bacteria_name, ".12hr.txt", sep=""), header=T)
    bacteria2 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_upregulated_DE_genes_in_", bacteria_name, ".36hr.txt", sep=""), header=T)
    bacteria3 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_upregulated_DE_genes_in_", bacteria_name, ".5.5d.txt", sep=""), header=T)
    bacteria4 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_downregulated_DE_genes_in_", bacteria_name, ".12hr.txt", sep=""), header=T)
    bacteria5 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_downregulated_DE_genes_in_", bacteria_name, ".36hr.txt", sep=""), header=T)
    bacteria6 = read.table(paste("/Users/JooHyun/Dropbox/RNA-Seq Project/5 - Analysis/Preliminary.analyses.Nov.2015/list_of_DEGs/list_of_downregulated_DE_genes_in_", bacteria_name, ".5.5d.txt", sep=""), header=T)
                                                                         
    colnames(bacteria1) = c("gene_id","gene_name","log2","FDR")
    colnames(bacteria2) = c("gene_id","gene_name","log2","FDR")
    colnames(bacteria3) = c("gene_id","gene_name","log2","FDR")
    colnames(bacteria4) = c("gene_id","gene_name","log2","FDR")
    colnames(bacteria5) = c("gene_id","gene_name","log2","FDR")
    colnames(bacteria6) = c("gene_id","gene_name","log2","FDR")
    
    bacteria = rbind(bacteria1, bacteria2, bacteria3, bacteria4, bacteria5, bacteria6); dim(bacteria)
    bacteria = unique(bacteria[,1]); length(bacteria)
    
    return(length(bacteria))
}
M.luteus_DE.counts = count_DEGs_in_each_condition_without_overlap("M.luteus"); M.luteus_DE.counts
E.coli_DE.counts = count_DEGs_in_each_condition_without_overlap("E.coli"); E.coli_DE.counts






