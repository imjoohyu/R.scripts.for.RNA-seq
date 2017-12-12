#piechart of core genes
#January 12th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/")

call_data = function(data_line){
    data <- read.table(data_line, sep ="\t", header=F)
    data = data[,2:5]
    colnames(data) = c("GO terms", "Number of genes", "Percent of gene hit against total number of genes", "Percent of gene hit against total number of process hits")
    data = data[order(data[,2],decreasing = T),]
    return(data)
}

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

#Call out a panel of 12 colors
color_panel = brewer.pal(12, "Set3")

#Upregulated genes
upregulated = call_data("panther_pie_chart_core_genes_7_plus_upregulated.txt")
upregulated_pie_chart = ggplot(upregulated, aes(x="", y=`Number of genes`, fill=`GO terms`)) + scale_fill_manual(values=color_panel) + geom_bar(width=1, stat="identity") + theme_minimal(base_size=30) + coord_polar("y", start=0) + theme(legend.position="bottom")
upregulated_pie_chart

#+ geom_label_repel(aes(label=`Percent of gene hit against total number of genes`))


#Downregulated genes
downregulated = call_data("panther_pie_chart_core_genes_7_plus_downregulated.txt")
downregulated_pie_chart = ggplot(downregulated, aes(x="", y=`Number of genes`, fill=`GO terms`)) + scale_fill_manual(values=c(color_panel[1],color_panel[2],color_panel[3],color_panel[4],color_panel[5],color_panel[7], color_panel[9],color_panel[10],color_panel[12])) + geom_bar(width=1, stat="identity") + coord_polar("y", start=0) + theme_minimal(base_size=30) + theme(legend.position="bottom")
downregulated_pie_chart
