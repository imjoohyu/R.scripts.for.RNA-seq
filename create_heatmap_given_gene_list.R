#Create a heatmap of any genes
#March 7th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

#Just a heatmap
gene_list = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others.txt", header=T, stringsAsFactors = F)
#gene_list = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others_12hr_only.txt", header=T, stringsAsFactors = F)
gene_list = gene_list[,1:2]
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #count data, not FC

expression_data_subset = expression_data[match(gene_list[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(4,7,10,13,16,19,22,25,26,27,28,29,32)] #12hr
rownames(expression_data_subset_with_name) = c(gene_list[,2])
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)


#pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEGs_up_in_some_conditions_and_down_in_others_heatmap.pdf", height=15, width=15) 
heatmap.2(as.matrix(expression_data_subset_with_name), Colv=FALSE, dendrogram="row", trace="none", symm=F, scale="row", col=my_palette, srtCol=45, key=T, keysize = 1, cexCol=1.2, cexRow=1, key.xlab="Expression")
#dev.off()


#Heatmap with a functional label
PC1_top_100_genes = read.table("PCA_using_DESeq2/genes.loaded.on.PC1.for.live.and.heatkilled.12hr_top100.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #count data, not FC

PC1_top_100_genes = PC1_top_100_genes[order(PC1_top_100_genes$gene_function),] #group by function
expression_data_subset = expression_data[match(PC1_top_100_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(4,7,10,13,16,19,22,25,26,27,28,29,32)]
rownames(expression_data_subset_with_name) = c(PC1_top_100_genes[,2])
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
other_colors = brewer.pal(12, "Set3")
table(PC1_top_100_genes[,3])
Label = c(other_colors[1],other_colors[2],rep(other_colors[3],31),rep(other_colors[4],2),rep(other_colors[5],6),rep(other_colors[6],2),rep(other_colors[7],3),rep(other_colors[8],2),rep(other_colors[9],52))

pdf(file="/Users/JooHyun/Desktop/t.pdf", height=5, width=5) 
heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", trace="none", symm=F, scale="row", col=my_palette, srtCol=45, key=T, keysize = 1, cexCol=0.76, cexRow=0.5, RowSideColors = Label, key.xlab="Expression")
dev.off()