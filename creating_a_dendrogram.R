#Create a dendrogram based on the original data
#February 28, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry

setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data_table <- read.table("edgeR_basic_comparison_FDR_converted_to_Y-N.txt",header=T) #11911
#library(gplots); library(RColorBrewer)
#my_palette = colorRampPalette(c("red","black","green"))(n = 299)

#Dendrogram based on 12hr sample >> This pattern goes well with the PCA pattern as expected (12/9/2016)
odd_index = seq(3,64,2)
data_table_FC_only = data_table[,c(odd_index)]; rownames(data_table_FC_only) = data_table[,1]
data_table_FC_only_12hr = data_table_FC_only[,c(1,4,7,10,13,16,19,22,23,24,25,26,29)] #12hr sample only, 11911 x 13
data_table_FC_only_12hr = t(data_table_FC_only_12hr) #13 x 11911

hc = hclust(dist(data_table_FC_only_12hr))
dd <- as.dendrogram(hc)
plot(hc, hang = -1); plot(rev(hc), hang = -1)

#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/dendrogram_of_conditions_out_of_11911.pdf", height=10, width=10) 
#heatmap.2(as.matrix(data_table_FC_only_12hr), Rowv=FALSE, density.info="none", dendrogram="column", trace="none", scale=c("column"), col=my_palette, srtCol=45)
#dev.off()
