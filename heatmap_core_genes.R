#plot a heatmap of core genes
#March 14th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

########################################
#Core upregulated genes

selected_core_genes = read.table("finding.core.genes/selected_core_genes_for_heatmap_after_re-annotation_040617.txt", header=T, stringsAsFactors = F)
#expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #count data, not FC
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data


expression_data_subset = expression_data[match(selected_core_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(9,15,21,27,33,39,45,47,49,51)] #for FC
#expression_data_subset_with_name = expression_data_subset[,c(7,10,13,16,19,22,25,26,27,28)] #for count
rownames(expression_data_subset_with_name) = c(selected_core_genes[,2])
colnames(expression_data_subset_with_name) = c("M.luteus", "E.coli", "S.marcescens Type", "E.faecalis", "P.rettgeri", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila")
expression_data_subset_with_name = expression_data_subset_with_name[,c(1:3,6,5,4,7:10)]

library(gplots); library(RColorBrewer); library("devtools") 
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
shades_for_bar = c("#08519c","#3182bd","#6baed6","#bdd7e7") #for the number of conditions
colors_for_bars = c("#a6611a", "#dfc27d", "#ffffbf", "#80cdc1", "#018571") #for functional groups
yes_no_bars = c("black") #for showing up only in live infections or not

#Create column color bars
selected_core_genes[,c(3:11)] = as.numeric(selected_core_genes[,c(3:11)])

#for the number of conditions
Numcondition = replace(selected_core_genes$N_sig, selected_core_genes$N_sig==10, shades_for_bar[1])
Numcondition = replace(Numcondition, Numcondition==9, shades_for_bar[2])
Numcondition = replace(Numcondition, Numcondition==8, shades_for_bar[3])
Numcondition = replace(Numcondition, Numcondition==7, shades_for_bar[4])

#for whether only in live infections or not
liveInfonly = replace(selected_core_genes$only_sig_in_live, selected_core_genes$only_sig_in_live==1, yes_no_bars[1])
liveInfonly = replace(liveInfonly, liveInfonly ==0, "white")

#for the functional groups
AntR = replace(selected_core_genes$Antimicrobial_Response, selected_core_genes$Antimicrobial_Response==1, colors_for_bars[1])
AntR = replace(AntR, AntR==0, "white")
StrR = replace(selected_core_genes$Stress_Response, selected_core_genes$Stress_Response==1, colors_for_bars[2])
StrR = replace(StrR, StrR==0, "white")
Metab = replace(selected_core_genes$Metabolism, selected_core_genes$Metabolism==1, colors_for_bars[3])
Metab = replace(Metab, Metab==0, "white")
Secr = replace(selected_core_genes$Secretion, selected_core_genes$Secretion==1, colors_for_bars[4])
Secr = replace(Secr, Secr==0, "white")
MetionR = replace(selected_core_genes$metal_ion_homeostasis, selected_core_genes$metal_ion_homeostasis==1, colors_for_bars[5])
MetionR = replace(MetionR, MetionR==0, "white")

#create a white bar for space (white regardless)
whiteBar = replace(selected_core_genes$TF, selected_core_genes$TF==1, "white")
whiteBar = replace(whiteBar, whiteBar==0, "white")

#put all of them together
rlab=t(cbind(Numcondition,whiteBar,liveInfonly, whiteBar,AntR, StrR, Metab, Secr, MetionR))
rownames(rlab)=c("NumC","WhiteBar","liveInf","WhiteBar", "AR", "SR", "M", "S", "I")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."


pdf(file="/Users/JooHyun/Desktop/t.pdf", height=14, width=14) 
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette, key=F, keysize = 3, cexCol=1.3, cexRow=0.75, RowSideColors=rlab, RowSideColorsSize=92) #gives a heatmap with an order based on PCA clustering
legend("topright",legend=c("Antimicrobial Response","Stress Response","Metabolism","Secretion","Metal Ion Homeostasis","","10","9","8","7","","Regulated only in live infection"),fill=c(colors_for_bars,"white",shades_for_bar,"white",yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()




########################################
#Core downregulated genes

selected_core_genes = read.table("finding.core.genes/selected_core_genes_for_heatmap_after_re-annotation_040617.txt", header=T, stringsAsFactors = F)
#expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #count data, not FC
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data


expression_data_subset = expression_data[match(selected_core_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(9,15,21,27,33,39,45,47,49,51)] #for FC
#expression_data_subset_with_name = expression_data_subset[,c(7,10,13,16,19,22,25,26,27,28)] #for count
rownames(expression_data_subset_with_name) = c(selected_core_genes[,2])
colnames(expression_data_subset_with_name) = c("M.luteus", "E.coli", "S.marcescens Type", "E.faecalis", "P.rettgeri", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila")
expression_data_subset_with_name = expression_data_subset_with_name[,c(1:3,6,5,4,7:10)]

library(gplots); library(RColorBrewer); library("devtools") 
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
shades_for_bar = c("#08519c","#3182bd","#6baed6","#bdd7e7") #for the number of conditions
colors_for_bars = c("#a6611a", "#dfc27d", "#ffffbf", "#80cdc1", "#018571") #for functional groups
yes_no_bars = c("black") #for showing up only in live infections or not

#Create column color bars
selected_core_genes[,c(3:11)] = as.numeric(selected_core_genes[,c(3:11)])

#for the number of conditions
Numcondition = replace(selected_core_genes$N_sig, selected_core_genes$N_sig==10, shades_for_bar[1])
Numcondition = replace(Numcondition, Numcondition==9, shades_for_bar[2])
Numcondition = replace(Numcondition, Numcondition==8, shades_for_bar[3])
Numcondition = replace(Numcondition, Numcondition==7, shades_for_bar[4])

#for whether only in live infections or not
liveInfonly = replace(selected_core_genes$only_sig_in_live, selected_core_genes$only_sig_in_live==1, yes_no_bars[1])
liveInfonly = replace(liveInfonly, liveInfonly ==0, "white")

#for the functional groups
AntR = replace(selected_core_genes$Antimicrobial_Response, selected_core_genes$Antimicrobial_Response==1, colors_for_bars[1])
AntR = replace(AntR, AntR==0, "white")
StrR = replace(selected_core_genes$Stress_Response, selected_core_genes$Stress_Response==1, colors_for_bars[2])
StrR = replace(StrR, StrR==0, "white")
Metab = replace(selected_core_genes$Metabolism, selected_core_genes$Metabolism==1, colors_for_bars[3])
Metab = replace(Metab, Metab==0, "white")
Secr = replace(selected_core_genes$Secretion, selected_core_genes$Secretion==1, colors_for_bars[4])
Secr = replace(Secr, Secr==0, "white")
MetionR = replace(selected_core_genes$metal_ion_homeostasis, selected_core_genes$metal_ion_homeostasis==1, colors_for_bars[5])
MetionR = replace(MetionR, MetionR==0, "white")

#create a white bar for space (white regardless)
whiteBar = replace(selected_core_genes$TF, selected_core_genes$TF==1, "white")
whiteBar = replace(whiteBar, whiteBar==0, "white")

#put all of them together
rlab=t(cbind(Numcondition,whiteBar,liveInfonly, whiteBar,AntR, StrR, Metab, Secr, MetionR))
rownames(rlab)=c("NumC","WhiteBar","liveInf","WhiteBar", "AR", "SR", "M", "S", "I")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."


pdf(file="/Users/JooHyun/Desktop/t.pdf", height=14, width=14) 
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette, key=F, keysize = 3, cexCol=1.3, cexRow=0.75, RowSideColors=rlab, RowSideColorsSize=92) #gives a heatmap with an order based on PCA clustering
legend("topright",legend=c("Antimicrobial Response","Stress Response","Metabolism","Secretion","Metal Ion Homeostasis","","10","9","8","7","","Regulated only in live infection"),fill=c(colors_for_bars,"white",shades_for_bar,"white",yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()
