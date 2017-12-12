###############################################################################
#Compare the list of DE genes to the DIRGs-385 genes
#August 24, 2016 (Updated on October 18th, 2016)
#Joo Hyun Im (ji72)
###############################################################################

#Remove previous entries
rm(list=ls(all=TRUE))

#Change:
#11/01/2016: There are actually 382 DIRGs, not 385 DIRGs. The DIRG_385_updated_simple.txt file's name got changed to DIRG_382_updated_2016.txt. The contents did not change.
#5/22/2017: I've changed the script so that the venn diagram can be drawn between DIRGs, core response, and response to heatkilled bacteria. I've also added the venn diagram with M.luteus DEGs and E.coli DEGs
#5/22/2017: I found a duplicate gene listed with different CG numbers in DIRG file. So I've updated it to DIRG_381_updated_2017.txt.


#Index:
#i) Compare the list of true core genes (8+) to the 385 DIRGs ==> Changed to 7+
#ii) Compare the list of all DEGs (1+ in at least one infection, including clean prick, infected, and heatkilled) to the 385 DIRGs
#iii) Compare the list of all DEGs (1+ in at least one live infection, excluding clean prick, infected, and heatkilled) to the 385 DIRGs = just live bacteria)



#================================================================================================
#i) Compare the list of true core genes (8+) to the 385 DIRGs ==> Changed to 7+
#full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes.up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
core.genes.up$type <- rep("Up", dim(core.genes.up)[1])
core.genes.down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
core.genes.down$type <- rep("Down", dim(core.genes.down)[1])
core.genes = rbind(core.genes.up, core.genes.down)

core.genes.names = matrix(core.genes[,1]) #252 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/DIRG_382_updated_2016.txt", header=T) #The original CG numbers from the Lemaitre document are converted into FBgn Flybase IDs (11/12/2015)
overlap = vector(); non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = dirg[which(dirg[,1] == core.genes.names[i,1]),]; #print(index)
    #cat("print i: ",i); cat(", index number: ", length(index)); print("\n")
    if (dim(index)[1] >= 1){ #if they overlap
        overlap = rbind(overlap,core.genes[i,])
    }
    else { #if they do not overlap
        non.overlap = rbind(non.overlap,core.genes[i,])
    }
}
dim(overlap) #shows the genes that overlap between my true core genes and DIRGs = 84 genes
dim(non.overlap) #shows the genes that do not overlap between my true core genes and DIRGs = 168 genes
overlap.2 = cbind(overlap,c(rep("Y",dim(overlap)[1]))); non.overlap.2 = cbind(non.overlap,c(rep("N",dim(non.overlap)[1])))
colnames(overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "regulation","DIRGs")
colnames(non.overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "regulation", "DIRGs")
core.genes.with.dirg = rbind(overlap.2, non.overlap.2); core.genes.with.dirg.sorted = core.genes.with.dirg[order(as.numeric(rownames(core.genes.with.dirg)), decreasing=F),]
write.table(core.genes.with.dirg.sorted, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7+_all_live_bacteria_compared_to_DIRGs.txt", quote=F, row.names=F, col.names=T)


#or simply
overlap_between_the_two = intersect(dirg[,1], core.genes[,1])


#========================================================================
#ii) Compare the list of all DEGs (1+ in at least one infection, including clean prick, infected, and heatkilled) to the 385 DIRGs
#full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
core.genes.names = matrix(core.genes[,1]) #2589 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/DIRG_382_updated_2016.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)

overlap = vector(); non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = dirg[which(dirg[,1] == core.genes.names[i,1]),]; #print(index)
    #cat("print i: ",i); cat(", index number: ", length(index)); print("\n")
    if (dim(index)[1] >= 1){ #if they overlap
        overlap = rbind(overlap,core.genes[i,])
    }
    else { #if they do not overlap
        non.overlap = rbind(non.overlap,core.genes[i,])
    }
}
dim(overlap) #shows the genes that overlap between my true core genes and DIRGs = 227 genes
dim(non.overlap) #shows the genes that overlap between my true core genes and DIRGs = 2362 genes
yes = matrix(c(rep("Y",dim(overlap)[1]))); no = matrix(c(rep("N",dim(non.overlap)[1]))); colnames(yes) = c("DIRGs"); colnames(no) = c("DIRGs")
overlap.2 = cbind(overlap[,c(1:2)], yes, overlap[,c(3:64)]); non.overlap.2 = cbind(non.overlap[,c(1:2)], no, non.overlap[,c(3:64)])
core.genes.with.dirg = rbind(overlap.2, non.overlap.2)
core.genes.with.dirg.sorted = core.genes.with.dirg[order(as.numeric(rownames(core.genes.with.dirg)), decreasing=F),]
write.table(core.genes.with.dirg.sorted, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEG-ups_in_all_conditions_compared_to_DIRGs.txt", quote=F, row.names=F, col.names=T)

#What are the 382-227 = 155 genes that are DIRGs but are not DEGs in our RNA-seq?
non.overlap.dirgs = vector()
for (i in 1:dim(dirg)[1]){
    index = core.genes.names[which(core.genes.names[,1] == as.character(dirg[i,1])),]
    if (length(index) == 0){ #if they do not overlap
        non.overlap.dirgs = rbind(non.overlap.dirgs, as.character(dirg[i,1]))
    }
}
dim(non.overlap.dirgs) #155

#Add the expression results to these non-overlapping genes
exp_results = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N.txt", header=T) #11911 x 64
non.overlap.dirgs.exp=c()
for (i in 1:dim(non.overlap.dirgs)[1]){
    index = exp_results[which(exp_results$gene_id == non.overlap.dirgs[i,1]),]
    non.overlap.dirgs.exp = rbind(non.overlap.dirgs.exp, index)
}

write.table(non.overlap.dirgs.exp, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/DIRGs_that_do_not_show_up_in_RNAseq.txt", quote=F, row.names=F, col.names=F)

#Add the expression results to these non-overlapping genes
exp_results_2 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", header=T) #17558 x 64
non.overlap.dirgs.exp2=c()
for (i in 1:dim(non.overlap.dirgs)[1]){
    index = exp_results_2[which(exp_results_2$gene.id == non.overlap.dirgs[i,1]),]
    non.overlap.dirgs.exp2 = rbind(non.overlap.dirgs.exp2,index)
}

write.table(non.overlap.dirgs.exp2, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/DIRGs_that_do_not_show_up_in_RNAseq_with_counts.txt", quote=F, row.names=F, col.names=T)




#========================================================================
#iii) Compare the list of all DEGs (1+ in at least one live infection, excluding clean prick, infected, and heatkilled) to the 385 DIRGs
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_num_of_conditions/core_genes_all_live_bacteria_upregulated.txt", header=T)
core.genes.names = matrix(core.genes[,1]) #1286 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_382_updated_2016.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
overlap = vector(); non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = dirg[which(dirg[,1] == core.genes.names[i,1]),]; #print(index)
    #cat("print i: ",i); cat(", index number: ", length(index)); print("\n")
    if (dim(index)[1] >= 1){ #if they overlap
        overlap = rbind(overlap,core.genes[i,])
    }
    else { #if they do not overlap
        non.overlap = rbind(non.overlap,core.genes[i,])
    }
}
dim(overlap) #shows the genes that overlap between my true core genes and DIRGs = 203 genes
dim(non.overlap) #shows the genes that overlap between my true core genes and DIRGs = 1080 genes
overlap.2 = cbind(overlap,c(rep("Y",dim(overlap)[1]))); non.overlap.2 = cbind(non.overlap,c(rep("N",dim(non.overlap)[1])))
colnames(overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "DIRGs")
colnames(non.overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "DIRGs")
core.genes.with.dirg = rbind(overlap.2, non.overlap.2); core.genes.with.dirg.sorted = core.genes.with.dirg[order(as.numeric(rownames(core.genes.with.dirg)), decreasing=F),]
write.table(core.genes.with.dirg.sorted, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_num_of_conditions/DEG-ups_in_all_live_conditions_compared_to_DIRGs.txt", quote=F, row.names=F, col.names=T)


all_genes_up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_num_of_conditions/core_genes_all_live_bacteria_upregulated.txt", header=T)
all_genes_down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_num_of_conditions/core_genes_all_live_bacteria_downregulated.txt", header=T)

all_genes_together = rbind(all_genes_up, all_genes_down) #2576
all_genes_together = unique(all_genes_together[,1])
overlap_with = intersect(all_genes_together, dirg[,1])




#Added on 5/22/2017

#Remove previous entries
rm(list=ls(all=TRUE))

#Comparison between DIRGs, core response, and response to heatkilled bacteria (added on 5/22/2017)
#core genes
core.genes.up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
core.genes.up$type <- rep("Up", dim(core.genes.up)[1])
core.genes.down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
core.genes.down$type <- rep("Down", dim(core.genes.down)[1])
core.genes = rbind(core.genes.up, core.genes.down) #252

#dirg
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DIRG/DIRG_381_updated_2017.txt", header=T) #The original CG numbers from the Lemaitre. 381 genes

#heatkilled
heatkilled = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.heatkilled.txt", header=F) #648


#get the overlaps
overlap_with_core_and_dirg = intersect(core.genes$gene_id, dirg$gene.id) #84
overlap_with_core_and_hk = intersect(core.genes$gene_id, heatkilled$V1) #142
overlap_with_dirg_and_hk = intersect(dirg$gene.id, heatkilled$V1) #77
total_overlap = intersect(overlap_with_core_and_dirg, heatkilled$V1) #59

core_and_dirg_minus_total_overlap = setdiff(overlap_with_core_and_dirg, total_overlap) #25
core_and_hk_minus_total_overlap = setdiff(overlap_with_core_and_hk, total_overlap) #83
dirg_and_hk_minus_total_overlap = setdiff(overlap_with_dirg_and_hk, total_overlap) #18

core.only = setdiff(core.genes$gene_id, overlap_with_core_and_dirg)
core.only = setdiff(core.only, core_and_hk_minus_total_overlap) #85
dirg.only = setdiff(dirg$gene.id, overlap_with_core_and_dirg)
dirg.only = setdiff(dirg.only, dirg_and_hk_minus_total_overlap) #279
hk.only = setdiff(heatkilled$V1, overlap_with_core_and_hk)
hk.only = setdiff(hk.only, dirg_and_hk_minus_total_overlap) #279

write.table(core.only, file="/Users/JooHyun/Desktop/overlap.txt", quote=F, row.names = F, col.names = F)

#the venn diagram with M.luteus DEGs and E.coli DEGs
#dirg
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DIRG/DIRG_381_updated_2017.txt", header=T) #The original CG numbers from the Lemaitre. 381 genes

#M. luteus DEGs
ml1 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_M.luteus.12hr.txt", header=T)
ml2 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_M.luteus.36hr.txt", header=T)
ml3 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_M.luteus.5.5d.txt", header=T)
ml4 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_M.luteus.12hr.txt", header=T)
ml5 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_M.luteus.36hr.txt", header=T)
ml6 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_M.luteus.5.5d.txt", header=T)
ml = c(as.character(ml1[,1]), as.character(ml2$gene_id), as.character(ml3$gene_id), as.character(ml4$gene_id), as.character(ml5$gene_id), as.character(ml6$gene_id))
ml_unique = unique(ml) #794

#E.coli DEGs
ec1 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_E.coli.12hr.txt", header=T)
ec2 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_E.coli.36hr.txt", header=T)
ec3 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_upregulated_DE_genes_in_E.coli.5.5d.txt", header=T)
ec4 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_E.coli.12hr.txt", header=T)
ec5 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_E.coli.36hr.txt", header=T)
ec6 = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/list_of_DEGs/list_of_downregulated_DE_genes_in_E.coli.5.5d.txt", header=T)
ec = c(as.character(ec1[,1]), as.character(ec2$gene_id), as.character(ec3$gene_id), as.character(ec4$gene_id), as.character(ec5$gene_id), as.character(ec6$gene_id))
ec_unique = unique(ec) #446

#get the overlaps
overlap_with_dirg_and_ml = intersect(dirg$gene.id, ml_unique) #124
overlap_with_dirg_and_ec = intersect(dirg$gene.id, ec_unique) #91
overlap_with_ml_and_ec = intersect(ml_unique, ec_unique) #297
total_overlap = intersect(overlap_with_dirg_and_ml, ec_unique) #80

dirg_and_ml_minus_total_overlap = setdiff(overlap_with_dirg_and_ml, total_overlap) #44
dirg_and_ec_minus_total_overlap = setdiff(overlap_with_dirg_and_ec, total_overlap) #11
ml_and_ec_minus_total_overlap = setdiff(overlap_with_ml_and_ec, total_overlap) #217

dirg.only = setdiff(dirg$gene.id, overlap_with_dirg_and_ml)
dirg.only = setdiff(dirg.only, dirg_and_ec_minus_total_overlap) #246
ml.only = setdiff(ml_unique, overlap_with_ml_and_ec)
ml.only = setdiff(ml.only, dirg_and_ml_minus_total_overlap) #453
ec.only = setdiff(ec_unique, overlap_with_ml_and_ec)
ec.only = setdiff(ec.only, dirg_and_ec_minus_total_overlap) #138

write.table(total_overlap, file="/Users/JooHyun/Desktop/overlap.txt", quote=F, row.names = F, col.names = F)



#Response to any live infection vs response to injury vs response to heat-killed bacteria (8/7/2017)
live_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T) #2423 genes
live_genes = data.frame(live_genes[,1]); colnames(live_genes) = c("V1") #2423 genes

heatkilled = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.heatkilled.txt", header=F) #648

injury_only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.cleanprick.txt", header=F) #114

live_genes = as.character(live_genes$V1)
heatkilled = as.character(heatkilled$V1)
injury_only = as.character(injury_only$V1)
length(intersect(live_genes, heatkilled))

int_btwn_live_injury = intersect(live_genes, injury_only)
int_btwn_live_hk = intersect(live_genes, heatkilled)
int_btwn_injury_hk = intersect(heatkilled, injury_only)

hk_1 = setdiff(heatkilled, int_btwn_live_hk)
hk_2 = setdiff(hk_1, int_btwn_injury_hk) #164 genes


#How close were these hk-specific genes in not making a cutoff for DE in live infections?
all_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T)
all_genes = all_genes[,c(1:52)] #remove heatkilled
hk_2 = data.frame(hk_2)
all_genes_selected = all_genes[which(all_genes$gene_id%in% hk_2$hk_2),]

#Count how many of the 164 genes had a p-value of <0.1 (close call) in live infections
all_genes_selected_edited = all_genes_selected
length.of.table <- as.numeric(dim(all_genes_selected)[1])
width.of.table <- as.numeric(dim(all_genes_selected)[2])
indicator <- NULL

for (k in 1:length.of.table){ #1, 2, 3, ... 11911
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        indicator <- isTRUE(all_genes_selected[k,m] > 0.1) #cutoff: FDR of 10%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) { #If the case is NOT significant,
            all_genes_selected_edited[k,m] = "N"
        }
        else {
            all_genes_selected_edited[k,m] = "Y"
        }
    }
}

write.table(all_genes_selected_edited, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/HK_specific_164_genes.txt", quote=F, col.names = T, row.names = F)


core.genes.up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
core.genes.up$type <- rep("Up", dim(core.genes.up)[1])
core.genes.down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
core.genes.down$type <- rep("Down", dim(core.genes.down)[1])
core_genes = rbind(core.genes.up, core.genes.down)
core_genes = as.character(core_genes$gene_id)

#Get the 493 genes that are HK specific when compared to core genes
int_btwn_core_injury = intersect(core_genes, injury_only)
int_btwn_core_hk = intersect(core_genes, heatkilled)
int_btwn_injury_hk = intersect(heatkilled, injury_only)

hk_1 = setdiff(heatkilled, int_btwn_core_hk)
hk_2 = setdiff(hk_1, int_btwn_injury_hk) #493 genes

all_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_all_genes_FC.txt", header=T)
all_genes = all_genes[,c(1:52)] #remove heatkilled
hk_2 = data.frame(hk_2)
all_genes_selected = all_genes[which(all_genes$gene_id%in% hk_2$hk_2),]

#Count how many of the 493 genes had a p-value of <0.1 (close call) in live infections
all_genes_selected_edited = all_genes_selected
length.of.table <- as.numeric(dim(all_genes_selected)[1])
width.of.table <- as.numeric(dim(all_genes_selected)[2])
indicator <- NULL

for (k in 1:length.of.table){ #1, 2, 3, ... 11911
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        indicator <- isTRUE(all_genes_selected[k,m] > 0.05) #cutoff: FDR of 5%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) { #If the case is NOT significant,
            all_genes_selected_edited[k,m] = "N"
        }
        else {
            all_genes_selected_edited[k,m] = "Y"
        }
    }
}

write.table(all_genes_selected_edited, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/HK_specific_493_genes_with_live_0.05.txt", quote=F, col.names = T, row.names = F)

