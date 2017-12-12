#Drawing a venn diagram that reflects the size of the samples
#Updated on 3/10/2017
#Joo Hyun Im (ji72)


#New script put together on 3/10/2017
library(VennDiagram); library(RColorBrewer)
draw_venn_diagram_two_groups = function(group1_value, group2_value, overlap_value, group1_name, group2_name){
    grid.newpage()
    draw.pairwise.venn(group1_value, group2_value, overlap_value, category = c(group1_name, group2_name), lty = rep("blank",2), fill = c("#ece7f2", "#2b8cbe"), scaled=T, cat.cex=1.3, cat.pos = c(0,0), cat.dist = rep(0.025, 2), fontface = "bold")
}
draw_venn_diagram_three_groups = function(group1_value, group2_value, group3_value, group12_value, group23_value, group13_value, overlap_value, group1_name, group2_name, group3_name){
    grid.newpage()
    draw.triple.venn(group1_value, group2_value, group3_value, group12_value, group23_value, group13_value, overlap_value, category = c(group1_name, group2_name, group3_name), lty = rep("blank",3), fill = c("#ece7f2", "#a6bddb","#2b8cbe"), scaled=T, cat.cex=2, fontface = "bold") 
}

#Gram-positive exclusive and Gram-negative exclusive
draw_venn_diagram_two_groups(662+1063, 851+1063, 1063, "Gram-positive exclusive","Gram-negative exclusive") 

#DIRGs and core genes
draw_venn_diagram_two_groups(168+84, 298+84, 84, "Core genes","DIRGs") 

#Core vs heatkilled vs sterile wound
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/core_heatkilled_stwound_venn_diagram.pdf", height=20, width=20) 
draw_venn_diagram_three_groups(105+5+78+64, 5+78+18+13, 64+78+13+493, 5+78, 78+13, 78+64, 78, "                     Core Response to 
                               Live Infection", "Response 
to Injury", "Response to MAMPs")
dev.off()

#Live vs heatkilled vs sterile wound
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/live_heatkilled_stwound_venn_diagram.pdf", height=20, width=20) 
draw_venn_diagram_three_groups(1918+22+90+393, 22+90+1+1, 393+90+1+164, 22+90, 90+1, 90+393, 90,"       Response to 
                         Live Infection","Response 
 to Injury", "Response to MAMPs")
dev.off()

#Core vs heatkilled vs sterile wound with simple text
pdf(file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/core_heatkilled_stwound_venn_diagram_v2.pdf", height=20, width=20) 
draw_venn_diagram_three_groups(105+5+78+64, 5+78+18+13, 64+78+13+493, 5+78, 78+13, 78+64, 78, "Core", "Injury", "MAMPs")
dev.off()


#Draw a venn diagram of DIRGs and core genes (old script prior to 3/10/2017)
# clean.prick = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs/list_of_upregulated_DE_genes_in_clean.prick_id_only.txt", header=T)
# clean.prick.uniq = as.matrix(unique(clean.prick)) #86
# core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 genes,
# #core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T) #2423 genes
# core.genes.names = matrix(core.genes[,1])  
# non.overlap = vector()
# for (i in 1:dim(core.genes.names)[1]){
#     index = clean.prick.uniq[which(clean.prick.uniq[,1] == core.genes.names[i,1]),]
#     if (length(index)<1){ #if there is no overlap
#         non.overlap = rbind(non.overlap,core.genes.names[i,1])
#     }
# }
# dim(non.overlap)
# 
# 
# library(VennDiagram)
# grid.newpage()
# # # for group 1, # for group 2, overlap
# 
# draw.pairwise.venn(386, 2589, 230, category = c("Group1", "Group 2"), 
#                    lty = rep("blank",2), fill = c("orange","skyblue"), alpha = rep(0.5, 2), cat.pos = c(0,0), 
#                    cat.dist = rep(0.025, 2), scaled=T) 
# #
# core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 genes
# core.genes.names = matrix(core.genes[,1]) #2589 genes
# dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_386_updated_simple.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
# overlap = vector(); non.overlap = vector()
# for (i in 1:dim(dirg)[1]){
#     index = core.genes.names[which(core.genes.names[,1] == dirg[i,1]),] #print(index)
#     #cat("print i: ",i); cat(", index number: ", length(index)); print("\n")
#     if (length(index) < 1){ #if they overlap
#         non.overlap = rbind(non.overlap, as.character(dirg[i,1]))
#     }
# }
# 
# write.table(overlap[,1:2], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/overlap.genes.from.DEGs.vs.DIRGs.txt", quote=F, row.names = F)
# write.table(cbind(non.overlap,output.mx), file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRGs.specific.genes.from.DEGs.vs.DIRGs.txt", quote=F, row.names = F)
