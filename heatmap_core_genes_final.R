#plot a heatmap of core genes
#May, 8th, 2017
#Joo Hyun Im (ji72)

#heatmap_core_genes.R was adjusted and saved here

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")


#Check which core genes are regulated only in live infections:
core = read.table("finding.core.genes/core_genes_for_heatmap/core_upgenes_list.txt", header=T, sep="\t") #166 genes
#core = read.table("finding.core.genes/core_genes_for_heatmap/core_dogenes_list.txt", header=T, sep="\t") #166 genes
live_only = read.table("specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.live.only.when.core_genes_were_used.txt", header=F) #105 genes

for(i in 1:dim(core)[1]){
    if ( as.character(core[i,1]) %in% live_only[,1] ){
        core[i,4] = 1
    }
    else{
        core[i,4] = 0
    }
}
write.table(core, file="finding.core.genes/core_genes_for_heatmap/up_sig_only_in_live_results.txt", quote=F, row.names = F, col.names = T)
#write.table(core, file="finding.core.genes/core_genes_for_heatmap/do_sig_only_in_live_results.txt", quote=F, row.names = F, col.names = T)




########################################
#Core upregulated genes
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
library(gplots); library(RColorBrewer); library("devtools") 

#Prepare data
selected_core_genes = read.table("finding.core.genes/core_genes_for_heatmap/core_upgenes_for_heatmap_function_classified_050817_simplified.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data
expression_data_subset = expression_data[match(selected_core_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(9,15,21,27,33,39,45,47,49,51)] #for FC
rownames(expression_data_subset_with_name) = c(selected_core_genes[,2])
colnames(expression_data_subset_with_name) = c("M.luteus", "E.coli", "S.marcescens Type", "E.faecalis", "P.rettgeri", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila")
expression_data_subset_with_name = expression_data_subset_with_name[,c(1:3,6,5,4,7:10)]

#Prepare colors
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)
#shades_for_bar = c("#525252","#969696", "#cccccc","#f7f7f7") #for the number of conditions (10, 9, 8, 7): grey
shades_for_bar = c("#525252","#969696", "#cccccc","#e2e2e2") #for the number of conditions (10, 9, 8, 7): grey, 7 is a bit darker here
colors_for_bars = c(brewer.pal(11,"Set3")[1:9], brewer.pal(11,"Set3")[11:12]) #for functional groups -- 11 functional groups
yes_no_bars = c("black") #for showing up only in live infections or not

#Create column color bars
selected_core_genes[,c(3:15)] = as.numeric(selected_core_genes[,c(3:15)])
selected_core_genes = selected_core_genes[,c(1:3,5:15)] #remove Signaling

#Correct the column names
colnames(selected_core_genes) = c("gene_id", "gene_name", "N.sig", "Antimicrobial response", "Stress response", "Metabolism", "Secretion", "Metal ion homeostasis", "Translation control", "Wound healing/tissue repair", "Cell redox/cell cycle control/cell growth", "Proteolysis", "Neuron-related","only_sig_in_live")


#for the number of conditions (10, 9, 8, 7)
Numcondition = replace(selected_core_genes$N.sig, selected_core_genes$N.sig==10, shades_for_bar[1])
Numcondition = replace(Numcondition, Numcondition==9, shades_for_bar[2])
Numcondition = replace(Numcondition, Numcondition==8, shades_for_bar[3])
Numcondition = replace(Numcondition, Numcondition==7, shades_for_bar[4])

#for whether only in live infections or not
liveInfonly = replace(selected_core_genes$only_sig_in_live, selected_core_genes$only_sig_in_live==1, yes_no_bars[1])
liveInfonly = replace(liveInfonly, liveInfonly ==0, "white")

#for the functional groups
C1 = replace(selected_core_genes[,4], selected_core_genes[,4]==1, colors_for_bars[1]); C1 = replace(C1, C1==0, "white")
C2 = replace(selected_core_genes[,5], selected_core_genes[,5]==1, colors_for_bars[2]); C2 = replace(C2, C2==0, "white")
C3 = replace(selected_core_genes[,6], selected_core_genes[,6]==1, colors_for_bars[3]); C3 = replace(C3, C3==0, "white")
C4 = replace(selected_core_genes[,7], selected_core_genes[,7]==1, colors_for_bars[4]); C4 = replace(C4, C4==0, "white")
C5 = replace(selected_core_genes[,8], selected_core_genes[,8]==1, colors_for_bars[5]); C5 = replace(C5, C5==0, "white")
C6 = replace(selected_core_genes[,9], selected_core_genes[,9]==1, colors_for_bars[6]); C6 = replace(C6, C6==0, "white")
C7 = replace(selected_core_genes[,10], selected_core_genes[,10]==1, colors_for_bars[7]); C7 = replace(C7, C7==0, "white")
C8 = replace(selected_core_genes[,11], selected_core_genes[,11]==1, colors_for_bars[8]); C8 = replace(C8, C8==0, "white")
C9 = replace(selected_core_genes[,12], selected_core_genes[,12]==1, colors_for_bars[9]); C9 = replace(C9, C9==0, "white")
C10 = replace(selected_core_genes[,13], selected_core_genes[,13]==1, colors_for_bars[10]); C10 = replace(C10, C10==0, "white")

#create a white bar for space (white regardless)
whiteBar = c(rep("white", each=length(C1)))

#put all of them together
rlab=t(cbind(Numcondition, whiteBar, liveInfonly, whiteBar, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10))

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."

#heatmap
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_for_heatmap/051117_final/core_upgenes_for_heatmap_17x35_ft1.9.pdf", height=35, width=17) 
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", key=F, col=my_palette, cexCol=1.3, cexRow = 1.9, RowSideColors=rlab, RowSideColorsSize=22, margin=c(3,9))
dev.off()

#legend
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", key=F, col=my_palette, cexCol=1.3, cexRow = 1.4, RowSideColors=rlab, RowSideColorsSize=22, margin=c(9,9))
#legend("topright",legend=c("Antimicrobial response","Stress response","Metabolism","Secretion","Metal ion homeostasis","Translation control", "Wound healing/tissue repair", "Cell redox/cell cycle control/cell growth", "Proteolysis", "Neuron-related","","10","9","8","7","","Regulated only in live infection"), fill=c(colors_for_bars,"white",shades_for_bar,"white",yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=1)
legend("topright",legend=c("10","9","8","7","","Regulated only in live infection"), fill=c(shades_for_bar,"white",yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

#colorkey
heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab="Expression Fold Change", cexCol=1.3, cexRow=0.75, margin=c(9,9), lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1))


########################################
#Core downregulated genes
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
library(gplots); library(RColorBrewer); library("devtools") 

#Prepare data
selected_core_genes = read.table("finding.core.genes/core_genes_for_heatmap/core_dogenes_for_heatmap_function_classified_050817_simplified.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data
expression_data_subset = expression_data[match(selected_core_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(9,15,21,27,33,39,45,47,49,51)] #for FC
rownames(expression_data_subset_with_name) = c(selected_core_genes[,2])
colnames(expression_data_subset_with_name) = c("M.luteus", "E.coli", "S.marcescens Type", "E.faecalis", "P.rettgeri", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila")
expression_data_subset_with_name = expression_data_subset_with_name[,c(1:3,6,5,4,7:10)]

#Prepare colors
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)
#shades_for_bar = c("#525252","#969696", "#cccccc","#f7f7f7") #for the number of conditions (10, 9, 8, 7): grey
shades_for_bar = c("#525252","#969696", "#cccccc","#e2e2e2") #for the number of conditions (10, 9, 8, 7): grey, 7 is a bit darker here
colors_for_bars = c(brewer.pal(11,"Set3")[11], brewer.pal(7,"Set2")[1], brewer.pal(11,"Set3")[3], brewer.pal(7,"Set2")[2:4], brewer.pal(11,"Set3")[9], brewer.pal(7,"Set2")[6:7])  #for functional groups -- 9 functional groups, colors matching with core upregulated
yes_no_bars = c("black") #for showing up only in live infections or not

#Create column color bars
selected_core_genes[,c(3:13)] = as.numeric(selected_core_genes[,c(3:13)])

#Correct the column names
colnames(selected_core_genes) = c("gene_id", "gene_name", "N.sig", "Neuron-related", "Response to toxin", "Metabolism", "Drug metabolism", "Reproduction", "Response to pheromone/olfactory behavior", "Proteolysis", "Oxidation reduction", "Carbohydrate binding", "only_sig_in_live")


#for the number of conditions (10, 9, 8, 7)
Numcondition = replace(selected_core_genes$N.sig, selected_core_genes$N.sig==10, shades_for_bar[1])
Numcondition = replace(Numcondition, Numcondition==9, shades_for_bar[2])
Numcondition = replace(Numcondition, Numcondition==8, shades_for_bar[3])
Numcondition = replace(Numcondition, Numcondition==7, shades_for_bar[4])

#for whether only in live infections or not
liveInfonly = replace(selected_core_genes$only_sig_in_live, selected_core_genes$only_sig_in_live==1, yes_no_bars[1])
liveInfonly = replace(liveInfonly, liveInfonly ==0, "white")

#for the functional groups
C1 = replace(selected_core_genes[,4], selected_core_genes[,4]==1, colors_for_bars[1]); C1 = replace(C1, C1==0, "white")
C2 = replace(selected_core_genes[,5], selected_core_genes[,5]==1, colors_for_bars[2]); C2 = replace(C2, C2==0, "white")
C3 = replace(selected_core_genes[,6], selected_core_genes[,6]==1, colors_for_bars[3]); C3 = replace(C3, C3==0, "white")
C4 = replace(selected_core_genes[,7], selected_core_genes[,7]==1, colors_for_bars[4]); C4 = replace(C4, C4==0, "white")
C5 = replace(selected_core_genes[,8], selected_core_genes[,8]==1, colors_for_bars[5]); C5 = replace(C5, C5==0, "white")
C6 = replace(selected_core_genes[,9], selected_core_genes[,9]==1, colors_for_bars[6]); C6 = replace(C6, C6==0, "white")
C7 = replace(selected_core_genes[,10], selected_core_genes[,10]==1, colors_for_bars[7]); C7 = replace(C7, C7==0, "white")
C8 = replace(selected_core_genes[,11], selected_core_genes[,11]==1, colors_for_bars[8]); C8 = replace(C8, C8==0, "white")
C9 = replace(selected_core_genes[,12], selected_core_genes[,12]==1, colors_for_bars[9]); C9 = replace(C9, C9==0, "white")

#create a white bar for space (white regardless)
whiteBar = c(rep("white", each=length(C1)))

#put all of them together
rlab=t(cbind(Numcondition, whiteBar, liveInfonly, whiteBar, C1, C2, C3, C4, C5, C6, C7, C8, C9))

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."

#heatmap
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_for_heatmap/051117_final/core_downgenes_for_heatmap_11x11_ft1.3.pdf", height=35, width=17) 
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", key=F, col=my_palette, cexCol=1.3, cexRow = 1.9, RowSideColors=rlab, RowSideColorsSize=22, margin=c(3,9))
dev.off()

#legend
heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", key=F, col=my_palette, cexCol=1.3, cexRow = 1.4, RowSideColors=rlab, RowSideColorsSize=22)
#legend("topright",legend=c("Neuron-related","Response to toxin", "Metabolism", "Drug metabolism", "Reproduction", "Response to pheromone/olfactory behavior", "Proteolysis", "Oxidation-reduction", "Carbohydrate binding","","","","","10","9","8","7","","Regulated only in live infection"), fill=c(colors_for_bars,"white",shades_for_bar,"white","white", "white", "white", yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=1)
legend("topright",legend=c("10","9","8","7","","Regulated only in live infection"), fill=c(shades_for_bar,"white",yes_no_bars), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

#colorkey
heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, Colv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab="Expression Fold Change", cexCol=1.3, cexRow=0.75, margin=c(9,9), lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1))
