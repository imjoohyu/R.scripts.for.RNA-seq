####################################################
#Find genes that load each PC in PCA
#September 7, 2016
#Joo Hyun Im (ji72)
####################################################

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")


####################################################
#1. Identify the genes in each PC
####################################################

#Using the CBSU machine
setwd("/workdir/ji72")
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103),]
design$treatment <- factor(design$treatment)

#From /Users/JooHyun/Dropbox/Cornell/Lab/R.scripts/deseq.with.time.R (Dec 2014/Jan 2015)
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable = read.table("mega_RNA-seq_count_of_all_samples.txt", header=T) #17558 x 105
design <- read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)
CountTable <- CountTable[,c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103)]
design <- design[c(8,11,14,17,20,23,26,27,28,29,30,33,43,46,49,52,55,58,61,62,63,64,65,68,78,81,84,87,90,93,96,97,98,99,100,103),]
design$treatment <- factor(design$treatment)

library("DESeq")
cds.all = newCountDataSet(CountTable, design)
#Normalization: estimate the effective library size.
cds.all = estimateSizeFactors(cds.all)
#Estimate variance without replicates
cdsBlind = estimateDispersions(cds.all, method="blind", sharingMode="fit-only") #if replicates are not available, estimate across conditions.

#PCA and then pull out the genes loading each PC
#The following method works similar to the one mentioned in https://www.biostars.org/p/13011/
library(genefilter)
vsd= varianceStabilizingTransformation(cdsBlind)
#vsd= varianceStabilizingTransformation(cds.all)
rv = rowVars(exprs(vsd))
select = order(rv, decreasing=TRUE)[seq_len(500)]
#select = order(rv, decreasing=TRUE) #getting everything, not just the top 500
pca = prcomp(t(exprs(vsd)[select,]))
fac = factor(apply(pData(vsd)[, c("treatment"), drop = FALSE], 1, paste, collapse = " : "))
pca.rotation <- as.matrix(pca$rotation[,1:2]) #list of genes with PC1 and PC2
write.table(pca.rotation, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt", quote=F, row.names = T, col.names = T)
#once I saved the list of genes loaded in each PC into a table,

PC1.top = data.frame(head(sort(pca.rotation[,1], decreasing = TRUE),20)); colnames(PC1.top) = c("PC1")
PC2.top = data.frame(head(sort(pca.rotation[,2], decreasing = TRUE),20)); colnames(PC2.top) = c("PC2")
write.table(PC1.top, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.top.PC1.for.live.and.heatkilled.12hr.txt", quote=F, row.names = T, col.names = T)
write.table(PC2.top, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.top.PC2.for.live.and.heatkilled.12hr.txt", quote=F, row.names = T, col.names = T)

PC1_top_from_17558_genes = sort(pca.rotation[,1], decreasing = TRUE)
#write.table(PC1_top_from_17558_genes, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr_all_17558_genes_PC1.txt", quote=F, row.names = T, col.names = T)


####################################################
#2. Are PC1 genes the Imd targets? (comparing the PC1 genes to the microarray data from De Gregorio et al 2002 -- see Evernote 3/29/17)
####################################################
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/")
PC1_genes = read.table("genes.loaded.on.PC1.for.live.and.heatkilled.12hr_top100_reannotated.txt", header=T, sep="\t") #100
#I first cleaned up the microarray data using sublime.
microarray_data = read.table("De_Gregorio_EMBO_2002_fully_processed_data_cleaned.txt", header=T, sep="\t") #13197 x 13

pick_Rel_targets = function(PC1_genes, microarray_data, cutoff){
    #Only pull the relevant microarray data
    microarray_data_subset = c()
    for (i in 1:dim(PC1_genes)[1]){
        microarray_data_subset = rbind(microarray_data_subset, microarray_data[grep(paste("\\b",PC1_genes$CG_name[i],"\\b", sep=""), microarray_data$CG_name),])
    }
    
    #Calculate Rel difference and WT difference
    microarray_data_subset$wt_diff = microarray_data_subset$wt90 - microarray_data_subset$wt0
    microarray_data_subset$rel_diff = microarray_data_subset$rel90 - microarray_data_subset$rel0
    
    for (j in 1:dim(microarray_data_subset)[1]){
        rel_over_wt =as.numeric(microarray_data_subset[j,15]/microarray_data_subset[j,14])
        if (rel_over_wt < cutoff){
            microarray_data_subset[j,16] = "Y"
        }
        else if (rel_over_wt <= 1){
            microarray_data_subset[j,16] = "N"
        }
        else if (rel_over_wt > 1){
            microarray_data_subset[j,16] = "Y"
        }
    }
    
    return(microarray_data_subset)
}
results = pick_Rel_targets(PC1_genes, microarray_data, 0.9)
#Saved the following and copied and pasted it into "genes.loaded.on.PC1.for.live.and.heatkilled.12hr_top100_reannotated.txt"
write.table(results[,c(1:3,6:7,14:16)], file="De_Gregorio_EMBO_2002_fully_processed_data_cleaned_with_wt_rel_comparison_cutoff_0.9.txt", quote=F, col.names=T, row.names=F)


####################################################
#3. Plot a heatmap of the PC1 genes
####################################################

#Plot a heatmap of expression from the top 100 genes from PC1 -- 12hr
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

PC1_top_100_genes = read.table("PCA_using_DESeq2/genes.loaded.on.PC1.for.live.and.heatkilled.12hr_top100_reannotated_w_rel_target_status.txt", header=T, stringsAsFactors = F, sep="\t")
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64
#expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #count data, not FC

library(gplots); library(RColorBrewer); library("devtools") 
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R") #when run, the following msg comes: "SHA-1 hash of file is.."
draw_a_heatmap = function(PC_genes, exp_data, order_by_what){
    cat("The genes in the heatmap are grouped by: ", order_by_what, "\n")
    
    #1) Get the data ready
    if (order_by_what == "function"){ #group the genes by function
        function_order=c("immune response", "metabolism", "proteolysis", "stress response", "wounding response", "neuronal activity", "metal ion transport", "cell redox homeostasis", "secretion", "NA")
        PC1_top_100_genes = PC1_top_100_genes[order(match(PC1_top_100_genes$gene_function, function_order)),] #group by function
    }
    else if (order_by_what == "Imd_target_status"){ #group the genes by imd target status
        function_order=c("Y", "N","NA")
        PC1_top_100_genes = PC1_top_100_genes[order(match(PC1_top_100_genes$rel_target_by_microarray, function_order)),] #group by target status
    }
    expression_data_subset = expression_data[match(PC1_top_100_genes[,1], expression_data[,1]),]
    expression_data_subset_with_name = expression_data_subset[,c(3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
    #expression_data_subset_with_name = expression_data_subset[,c(4,7,10,13,16,19,22,25,26,27,28,29,32)] #for count
    rownames(expression_data_subset_with_name) = PC1_top_100_genes$gene_name
    colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")
    
    #2) Get the colors ready
    my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
    other_colors = brewer.pal(12, "Set3")

    #3) Plot a histogram
    if (order_by_what == "function"){ #group the genes by function
        Label = c(rep("#a6611a",38),rep("#ffffbf",12),rep(other_colors[3],5),rep("#dfc27d",7),rep(other_colors[4],5),rep(other_colors[5],5),rep("#018571",3), other_colors[6], "#80cdc1", rep(other_colors[7],23)) #to be consistent with the core gene function colors #by function
        
        pdf(file="PCA_using_DESeq2/heatmap_PC1_genes_with_functional_annotation_RPlot.pdf", height=20, width=20) 
        heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="column",trace="none", symm=F, scale="none", col=my_palette, srtCol=45, key=T, keysize = 1, cexCol=1.5, cexRow=1, RowSideColors = Label, key.xlab="Expression") #gives a heatmap with an order based on PCA clustering
        dev.off()
    }
    
    else if (order_by_what == "Imd_target_status"){ #group the genes by imd target status
        PC1_top_100_genes$rel_target_by_microarray = as.character(PC1_top_100_genes$rel_target_by_microarray)
        PC1_top_100_genes$rel_target_by_microarray = replace(PC1_top_100_genes$rel_target_by_microarray, is.na(PC1_top_100_genes$rel_target_by_microarray), "NA")
        PC1_top_100_genes$gene_functiob = as.character(PC1_top_100_genes$gene_function)
        PC1_top_100_genes$gene_function = replace(PC1_top_100_genes$gene_function, is.na(PC1_top_100_genes$gene_function), "NA")
        
        Rel_target = replace(PC1_top_100_genes$rel_target_by_microarray, PC1_top_100_genes$rel_target_by_microarray=="Y", "black")
        Rel_target = replace(Rel_target, Rel_target=="N", "white")
        Rel_target = replace(Rel_target, Rel_target=="NA", "grey")
        
        gene_function_list = replace(PC1_top_100_genes$gene_function, PC1_top_100_genes$gene_function == "immune response", "#a6611a")
        gene_function_list = replace(gene_function_list, gene_function_list == "metabolism", "#ffffbf")
        gene_function_list = replace(gene_function_list, gene_function_list == "proteolysis", other_colors[3])
        gene_function_list = replace(gene_function_list, gene_function_list == "stress response", "#dfc27d")
        gene_function_list = replace(gene_function_list, gene_function_list == "wounding response", other_colors[4])
        gene_function_list = replace(gene_function_list, gene_function_list == "neuronal activity", other_colors[5])
        gene_function_list = replace(gene_function_list, gene_function_list == "metal ion transport", "#018571")
        gene_function_list = replace(gene_function_list, gene_function_list == "cell redox homeostasis", other_colors[6])
        gene_function_list = replace(gene_function_list, gene_function_list == "secretion", "#80cdc1")
        gene_function_list = replace(gene_function_list, gene_function_list == "NA", other_colors[7])
        
        rlab=t(cbind(Rel_target, gene_function_list))
        rownames(rlab)=c("Imd target","Gene function")
        
        pdf(file="PCA_using_DESeq2/heatmap_PC1_genes_with_functional_annotation_and_Imd_target_status_RPlot.pdf", height=16, width=20) 
        mai=c(2,2,1,3)
        plot = heatmap.3(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="column", trace="none", symm=F, scale="none", col=my_palette, key=F, keysize = 1, cexCol=1.2, cexRow=1, RowSideColors=rlab, RowSideColorsSize=92) #gives a heatmap with an order based on PCA clustering
        legend("topleft",legend=c("Imd target", "Non Imd target", "No imformation available", "Immune response", "Metabolism", "Proteolysis", "Stress response", "Wounding response", "Neuronal activity", "Metal ion transport", "Cell redox homeostasis", "Secretion", "Unknown"), fill=c("black","white","grey","#a6611a", "#ffffbf", other_colors[3],"#dfc27d", other_colors[4],other_colors[5], "#018571",other_colors[6],"#80cdc1", other_colors[7]), border=FALSE, bty="n", y.intersp = 0.7, cex=1.2)
        dev.off()
    }
}


#draw_a_heatmap(PC1_top_100_genes, expression_data, "function")
draw_a_heatmap(PC1_top_100_genes, expression_data, "Imd_target_status")



####################################################
#4. Plot a heatmap of the PC2 genes
####################################################
PC2_top_100_genes = read.table("PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt", header=T, stringsAsFactors = F)
expression_data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T) #FC, not count data, 11911 x 64

PC2_top_100_genes = PC2_top_100_genes[order(PC2_top_100_genes$PC2, decreasing = T),] #order by contribution



expression_data_subset = expression_data[match(rownames(PC2_top_100_genes), expression_data[,1]),]
expression_data_subset_with_name = expression_data_subset[,c(1,2,3,9,15,21,27,33,39,45,47,49,51,53,59)] #for FC
expression_data_subset_with_name = expression_data_subset_with_name[1:100,]
rownames(expression_data_subset_with_name) = expression_data_subset_with_name[,2]
expression_data_subset_with_name = expression_data_subset_with_name[,c(3:15)]
colnames(expression_data_subset_with_name) = c("Sterile Wound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live", "Ecc15", "S.aureus", "P.sneebia", "S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled", "P.rettgeri heatkilled")

library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
heatmap.2(as.matrix(expression_data_subset_with_name), Rowv=FALSE, density.info="none", dendrogram="column",trace="none", symm=F, scale="none", col=my_palette, srtCol=45, key=T, keysize = 1, cexCol=1.5, cexRow=1, key.xlab="Expression") #gives a heatmap with an order based on PCA clustering







#========What didn't work: the loadings were done with sample IDs not with gene ids========== (old work, updated as above on 10/4/2016)
# #DESeq2: I have to run DESeq2 because the PCA function below takes DESeq2 function 'assay()' which is part of the SummarizedExperiment package that DESeq2 automatically loads.
# library("DESeq2")
# design.treatment.only = data.frame(design[,-2:-3]) #make the design file succinct with just the treatment information
# rownames(design.treatment.only) = rownames(design); colnames(design.treatment.only) = c("treatment")
# design.treatment.only$treatment <- factor(design.treatment.only$treatment)
# attach(design.treatment.only)
# #construct the dataset from matrix
# cds=DESeqDataSetFromMatrix(countData=CountTable, colData = design.treatment.only, design = ~ treatment)
# cds = estimateSizeFactors(cds)
# cdsB = estimateDispersions(cds)
# vsd = varianceStabilizingTransformation(cdsB)
# 
# #Pull out appropriate values and print PCA plot
# #library(RColorBrewer)
# library(genefilter)
# library(lattice)
# #library(SummarizedExperiment)
# plotPCAWithSampleNames = function(x, intgroup=c("treatment"), ntop=500) {
#     # pca
#     rv = rowVars(assay(x))
#     select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
#     pca = prcomp(t(assay(x)[select,]))
#     
#     # proportion of variance
#     variance = pca$sdev^2 / sum(pca$sdev^2)
#     variance = round(variance, 3) * 100
#     
#     # sample names
#     names = colnames(x)
#     
#     # factor of groups
#     fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
#     
#     #     # colors
#     #     if( nlevels(fac) >= 10 )
#     #         colors = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"), brewer.pal(12,"Set3"), brewer.pal(7,"Dark2"))
#     #     
#     #     # plot
#     #     xyplot(
#     #         PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
#     #         aspect = "fill",
#     #         col = colors,
#     #         xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
#     #         ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
#     #         panel = function(x, y, ...) {
#     #             panel.xyplot(x, y, ...);
#     #             ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
#     #         },
#     #         main = draw.key(
#     #             key = list(
#     #                 rect = list(col = colors),
#     #                 text = list(levels(fac)),
#     #                 rep = FALSE
#     #             )
#     #         )
#     #     )
#     #     
#     cat("PC1 %: ", variance[1], ", PC2 %: ",variance[2], ", PC3 %: ",variance[3])
#     pca.data = as.data.frame(pca$x) #Pull out the genes loading each PC
#     write.table(pca.data[,c(1:3)], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt", quote=F, col.names=T,row.names = T)
# } #cut at the top 500 genes
# plotPCAWithSampleNames(vsd, intgroup=c("treatment"), ntop=500)
# # pdf(file="/home/ji72/RNAseq/totalRNAseq_Aug2014_hiseq/bam/DESeq_GLM_PCA.pdf", height=5, width=5) 
# # plotPCA(vsd, intgroup=c("species","time"))
# # dev.off()
# 
# 
# #once I saved the list of genes loaded in each PC into a table,
# pca.loaded <- read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/genes.loaded.on.PC1.and.PC2.for.live.and.heatkilled.12hr.txt")
# pca.loaded.sorted = pca.loaded[order(pca.loaded$PC1, decreasing = TRUE),]




