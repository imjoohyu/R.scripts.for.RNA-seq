#Creating a heatmap for PC1 and PC2 genes
#May 8th, 2017
#Joo Hyun Im (ji72)

####################################################
#1. Plot a heatmap of the PC1 genes
#Plot a heatmap of expression from the top 100 genes from PC1 -- 12hr
####################################################

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

PC1_top_100_genes = read.table("PCA_using_DESeq2/genes.loaded.on.PC1.for.live.and.heatkilled.12hr_top100_reannotated_w_rel_target_status.txt", header=T, stringsAsFactors = F, sep="\t")
PC1_top_100_genes = PC1_top_100_genes[-which(PC1_top_100_genes$gene_id == "FBgn0000280"),] #removing Cec-psi1 and adding CG30026.
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64

library(gplots); library(RColorBrewer); library("devtools") 
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."

#1) get the data ready
function_order=c("Y", "N","NA")
PC1_top_100_genes = PC1_top_100_genes[order(match(PC1_top_100_genes$rel_target_by_microarray_042717, function_order)),] #group by target status
expression_data_subset = expression_data[match(PC1_top_100_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
rownames(expression_data_subset_with_name) = PC1_top_100_genes$gene_name
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

#2) Get the colors ready for the color key
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)

#3) Plot a histogram
PC1_top_100_genes$rel_target_by_microarray_042717 = as.character(PC1_top_100_genes$rel_target_by_microarray_042717)
PC1_top_100_genes$rel_target_by_microarray_042717 = replace(PC1_top_100_genes$rel_target_by_microarray_042717, is.na(PC1_top_100_genes$rel_target_by_microarray_042717), "NA")

Rel_target = replace(PC1_top_100_genes$rel_target_by_microarray_042717, PC1_top_100_genes$rel_target_by_microarray_042717=="Y", "black")
Rel_target = replace(Rel_target, Rel_target=="N", "antiquewhite2") #used to be "white"
Rel_target = replace(Rel_target, Rel_target=="NA", "grey")

rlab=t(Rel_target)
rownames(rlab)=c("Imd target")

heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab= "Expression Fold Change", cexCol=1.2, cexRow=1.2, RowSideColors=rlab, RowSideColorsSize=10, margin = c(10,10)) #gives a heatmap with an order based on PCA clustering
legend("topleft",legend=c("Required by Imd", "Not required by Imd", "No information available"), fill=c("black","antiquewhite2","grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=1.2)
heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab = "Expression Fold Change", margin = c(10,10)) #just get the color key


####################################################
#2. Plot a heatmap of the PC2 genes
#Plot a heatmap of expression from the top 100 genes from PC2 -- 12hr
####################################################

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

PC2_top_100_genes = read.table("PCA_using_DESeq2/PC2_top100_genes_gap_filled_more_added.txt", header=T, stringsAsFactors = F) #already sorted by the amount of contribution
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0053462"),] #removing CG33462
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0000281"),] #removing Cec2
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0262722"),] #removing CR43166
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0262983"),] #removing CG43291
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0050272"),] #removing MFS1
PC2_top_100_genes = PC2_top_100_genes[-which(PC2_top_100_genes$gene_id == "FBgn0038214"),] #removing CG9616
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64

expression_data_subset = expression_data[match(PC2_top_100_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(1,2,3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
expression_data_subset_with_name = expression_data_subset_with_name[1:100,]
expression_data_subset_with_name$gene_name = as.character(droplevels.factor(expression_data_subset_with_name$gene_name))
rownames(expression_data_subset_with_name) = c(PC2_top_100_genes[,2])

expression_data_subset_with_name = expression_data_subset_with_name[,c(3:15)]
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

library(gplots); library(RColorBrewer); library(RColorBrewer)
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)

heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="column",trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab="Expression Fold Change", cexCol=1.2, cexRow=1, margin=c(10,10)) #gives a heatmap with an order based on PCA clustering




####################################################
#3. Plot a heatmap of the PC1 genes with Toll and Imd regulation noted
#Plot a heatmap of expression from the top 100 genes from PC1 -- 12hr
####################################################

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

PC1_top_100_genes = read.table("PCA_using_DESeq2/final_PC1_PC2_heatmaps/final_051017/PC1_top100_with_microarray_established.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64

library(gplots); library(RColorBrewer); library("devtools") 
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."

#1) get the data ready
#function_order=c("Y", "N","NA")
#PC1_top_100_genes = PC1_top_100_genes[order(match(PC1_top_100_genes$Toll_spz_cutoff, function_order)),] #group by target status
expression_data_subset = expression_data[match(PC1_top_100_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
rownames(expression_data_subset_with_name) = PC1_top_100_genes$gene_name
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

#2) Get the colors ready for the color key
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)

#3) Plot a histogram
PC1_top_100_genes$Toll_spz_cutoff =  as.character(PC1_top_100_genes$Toll_spz_cutoff)
PC1_top_100_genes$Toll_spz_cutoff = replace(PC1_top_100_genes$Toll_spz_cutoff, is.na(PC1_top_100_genes$Toll_spz_cutoff), "NA")
PC1_top_100_genes$Imd_rel_cutoff = as.character(PC1_top_100_genes$Imd_rel_cutoff)
PC1_top_100_genes$Imd_rel_cutoff = replace(PC1_top_100_genes$Imd_rel_cutoff, is.na(PC1_top_100_genes$Imd_rel_cutoff), "NA")

Toll_target = replace(PC1_top_100_genes$Toll_spz_cutoff, PC1_top_100_genes$Toll_spz_cutoff=="Y", "black")
Toll_target = replace(Toll_target, Toll_target=="N", "antiquewhite2") #used to be "white"
Toll_target = replace(Toll_target, Toll_target=="NA", "grey")
Imd_target = replace(PC1_top_100_genes$Imd_rel_cutoff, PC1_top_100_genes$Imd_rel_cutoff=="Y", "black")
Imd_target = replace(Imd_target, Imd_target=="N", "antiquewhite2") #used to be "white"
Imd_target = replace(Imd_target, Imd_target=="NA", "grey")

rlab=t(cbind(Toll_target, Imd_target))
rownames(rlab)=c("Regulated by Toll", "Regulated by Imd")

heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=F, keysize = 1, key.xlab= "Expression Fold Change", cexCol=1.2, cexRow=1.2, RowSideColors=rlab, RowSideColorsSize=10, margin = c(10,10)) #gives a heatmap with an order based on PCA clustering

heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab = "Expression Fold Change", margin = c(10,10)) #just get the color key
legend("topleft",legend=c("Required", "Not required", "No information available"), fill=c("black","antiquewhite2","grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=1.2) #legend


#For PC2 (8/3/2017)
m(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

PC2_top_100_genes = read.table("PCA_using_DESeq2/final_PC1_PC2_heatmaps/final_051017/PC2_top100_with_microarray_established.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64

library(gplots); library(RColorBrewer); library("devtools") 
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."

#1) get the data ready
expression_data_subset = expression_data[match(PC2_top_100_genes[,1], expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
rownames(expression_data_subset_with_name) = PC2_top_100_genes$gene_name
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

#2) Get the colors ready for the color key
my_palette = colorRampPalette(c("#4575b4","white","#d73027"))(n = 299) #blue, white, and red (final)

#3) Plot a histogram
PC2_top_100_genes$Toll_spz_cutoff =  as.character(PC2_top_100_genes$Toll_spz_cutoff)
PC2_top_100_genes$Toll_spz_cutoff = replace(PC2_top_100_genes$Toll_spz_cutoff, is.na(PC2_top_100_genes$Toll_spz_cutoff), "NA")
PC2_top_100_genes$Imd_rel_cutoff = as.character(PC2_top_100_genes$Imd_rel_cutoff)
PC2_top_100_genes$Imd_rel_cutoff = replace(PC2_top_100_genes$Imd_rel_cutoff, is.na(PC2_top_100_genes$Imd_rel_cutoff), "NA")

Toll_target = replace(PC2_top_100_genes$Toll_spz_cutoff, PC2_top_100_genes$Toll_spz_cutoff=="Y", "black")
Toll_target = replace(Toll_target, Toll_target=="N", "antiquewhite2") #used to be "white"
Toll_target = replace(Toll_target, Toll_target=="NA", "grey")
Imd_target = replace(PC2_top_100_genes$Imd_rel_cutoff, PC2_top_100_genes$Imd_rel_cutoff=="Y", "black")
Imd_target = replace(Imd_target, Imd_target=="N", "antiquewhite2") #used to be "white"
Imd_target = replace(Imd_target, Imd_target=="NA", "grey")

rlab=t(cbind(Toll_target, Imd_target))
rownames(rlab)=c("Regulated by Toll", "Regulated by Imd")

heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=F, keysize = 1, key.xlab= "Expression Fold Change", cexCol=1.2, cexRow=1.2, RowSideColors=rlab, RowSideColorsSize=10, margin = c(10,10)) #gives a heatmap with an order based on PCA clustering

heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="none", trace="none", symm=F, scale="none", col=my_palette, key=T, keysize = 1, key.xlab = "Expression Fold Change", margin = c(10,10)) #just get the color key
legend("topleft",legend=c("Required", "Not required", "No information available"), fill=c("black","antiquewhite2","grey"), border=FALSE, bty="n", y.intersp = 0.7, cex=1.2) #legend


