#Plot genes based on expression level and chromosome positions
#January 13th, 2017
#Joo Hyun Im (ji72)


rm(list=ls(all=TRUE)) #delete any previous entry

setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")
data_table = read.table("edgeR_basic_comparison_all_genes_FC.txt",header=T)
coordinates = read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt", header=T)

#1. Match the chr positions to the genes
match_the_coor = function(data_table){
    coor_table = matrix(NA, nrow=dim(data_table)[1], ncol=3)
    for (i in 1:dim(data_table)[1]){
        gene_id = as.character(data_table[i,1])
        chr = coordinates[which(coordinates$gene_id == gene_id),1]
        coor_1 = coordinates[which(coordinates$gene_id == gene_id),4]
        coor_2 = coordinates[which(coordinates$gene_id == gene_id),5]
        coor_table[i,1] = as.character(chr)
        coor_table[i,2] = as.numeric(coor_1)
        coor_table[i,3] = as.numeric(coor_2)
    }
    return(coor_table)
}
coordinates_matched = match_the_coor(data_table)

#2. Plot the genes
library(ggplot2); library(grid); library(gridExtra); library(grepel)
data_table_coor_matched = cbind(data_table[,c(1,2,13,14)], coordinates_matched)
colnames(data_table_coor_matched) = c("gene_id","gene_name","log2FC","FDR","chr", "coor_1", "coor_2")

chr_list = c("2L","2R","3L","3R","4","X")
grid.newpage(); pushViewport(viewport(layout = grid.layout(6, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
for(j in 1:length(chr_list)){
    data_table_coor_matched_chr = data_table_coor_matched[which(data_table_coor_matched$chr == chr_list[j]),]
    data_table_coor_matched_chr$coor_1 = as.numeric(levels(data_table_coor_matched_chr$coor_1))[data_table_coor_matched_chr$coor_1]
    data_table_coor_matched_chr$coor_2 = as.numeric(levels(data_table_coor_matched_chr$coor_2))[data_table_coor_matched_chr$coor_2]
    data_table_coor_matched_chr_sorted = data_table_coor_matched_chr[order(data_table_coor_matched_chr$coor_1, decreasing = F),]
    for(i in 1:dim(data_table_coor_matched_chr_sorted)[1]){
        if (data_table_coor_matched_chr_sorted[i,4] < 0.05){
            data_table_coor_matched_chr_sorted[i,4] = "TRUE"
        }
        else{
            data_table_coor_matched_chr_sorted[i,4] = "FALSE"
        }
    }

    plot = ggplot(data_table_coor_matched_chr_sorted, aes(x=coor_1, y=log2FC, fill=FDR)) + geom_point(aes(color=FDR), alpha=3/4) + theme_minimal(base_size=15) + scale_colour_manual(values = c("lightgrey","purple")) + labs(title=paste("Chr",chr_list[j],sep=""))
    print(plot, vp = vplayout(j, 1))
#     if (j<3){
#         print(plot, vp = vplayout(1, j))
#     }
#     else if (j==3 | j==4){
#         print(plot, vp = vplayout(2, j-2))
#     }
#     else{
#         print(plot, vp = vplayout(3, j-4))
#     }
}

