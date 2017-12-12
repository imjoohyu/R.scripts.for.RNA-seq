clean.prick = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs/list_of_upregulated_DE_genes_in_clean.prick_id_only.txt", header=T)
clean.prick.uniq = as.matrix(unique(clean.prick)) #86
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 genes,
#core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T) #2423 genes
core.genes.names = matrix(core.genes[,1])  
non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = clean.prick.uniq[which(clean.prick.uniq[,1] == core.genes.names[i,1]),]
    if (length(index)<1){ #if there is no overlap
        non.overlap = rbind(non.overlap,core.genes.names[i,1])
    }
}
dim(non.overlap)


library(VennDiagram)
grid.newpage()

# # for group 1, # for group 2, overlap

draw.pairwise.venn(386, 2589, 230, category = c("Group1", "Group 2"), 
                   lty = rep("blank",2), fill = c("orange","skyblue"), alpha = rep(0.5, 2), cat.pos = c(0,0), 
                   cat.dist = rep(0.025, 2), scaled=T) 


#
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 genes
core.genes.names = matrix(core.genes[,1]) #2589 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_386_updated_simple.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
overlap = vector(); non.overlap = vector()
for (i in 1:dim(dirg)[1]){
    index = core.genes.names[which(core.genes.names[,1] == dirg[i,1]),] #print(index)
    #cat("print i: ",i); cat(", index number: ", length(index)); print("\n")
    if (length(index) < 1){ #if they overlap
        non.overlap = rbind(non.overlap, as.character(dirg[i,1]))
    }
}

write.table(overlap[,1:2], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/overlap.genes.from.DEGs.vs.DIRGs.txt", quote=F, row.names = F)
write.table(cbind(non.overlap,output.mx), file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRGs.specific.genes.from.DEGs.vs.DIRGs.txt", quote=F, row.names = F)
