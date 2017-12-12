#Pulling out DE genes for a given condition in Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_FCgene_ID-name_FDR_converted_to_Y-N_at_least_one_sig.txt 
#October 12th, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")
data.table <- read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt",header=T)
#data.table = read.table("Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_FCgene_ID-name_with_pval_at_least_one_sig.txt", header=T)

list.of.name.of.conditions <- list("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d",
                                   "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                   "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
                                   "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")

#1. Pull out DE genes for a given condition
count=1
for (n in seq(4,64,2)){
    sig.genes.in.this.condition = subset(data.table, grepl("Y",data.table[,n]) ); print(dim(sig.genes.in.this.condition)[1])
    subset = sig.genes.in.this.condition[,c(1:2,n-1,n)]
    subset.sorted = subset[order(subset[,3], decreasing=T),]
    write.table(subset.sorted[which(subset.sorted[,3] > 0),], file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_upregulated_DE_genes_in_",list.of.name.of.conditions[count],".txt",sep=""), 
                quote=F, col.names= T, row.names=F)
    write.table(subset.sorted[which(subset.sorted[,3] < 0),], file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_downregulated_DE_genes_in_",list.of.name.of.conditions[count],".txt",sep=""), 
                quote=F, col.names= T, row.names=F)
    #write.table(subset, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/GO_analysis_by_DAVID_Oct_2015/list_of_DE_genes_in_",list.of.name.of.conditions[count],".txt",sep=""),
    #           quote=F, col.names= T, row.names=F)
    count = count +1
}

#2. DAVID on https://david.ncifcrf.gov/ -- Save the output as one of the names in GO.list below

#3. Process the DAVID outputs
####### Update: no longer in use as of 2/22/2017 (Wed)
# setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA")
# 
# direction = c("up","down") #For downregulated genes, change the number for direction to direction[2] and get rid of "GO_clean.prick.5.5d",GO_clean.prick.36hr","GO_E.coli.36hr"   from the list
# GO.list = c("GO_clean.prick.12hr", "GO_clean.prick.36hr", "GO_clean.prick.5.5d", 
#             "GO_M.luteus.12hr", "GO_M.luteus.36hr", "GO_M.luteus.5.5d", "GO_E.coli.12hr", "GO_E.coli.36hr", 
#             "GO_E.coli.5.5d","GO_S.mar.type.12hr", "GO_S.mar.type.36hr", "GO_S.mar.type.5.5d", "GO_E.fae.live.12hr", "GO_E.fae.live.36hr", "GO_E.fae.live.5.5d", "GO_P.rett.live.12hr", "GO_P.rett.live.36hr",
# "GO_P.rett.live.5.5d", "GO_Ecc15.12hr", "GO_Ecc15.36hr", "GO_Ecc15.5.5d", "GO_S.aureus.12hr", "GO_P.sneebia.12hr", "GO_S.mar.Db11.12hr", "GO_P.ento.12hr",
# "GO_E.fae.heatkilled.12hr", "GO_E.fae.heatkilled.36hr", "GO_E.fae.heatkilled.5.5d", "GO_P.rett.heatkilled.12hr", "GO_P.rett.heatkilled.36hr", "GO_P.rett.heatkilled.5.5d")
# 
# list.of.name.of.conditions <- list("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", 
#                                    "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", 
#                                    "E.coli.5.5d", "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
#                                    "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
#                                    "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
# 
# j = 1; #j = 2
# for (m in 1:length(GO.list)){
#     cat("Reading in the following file: ", paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/", GO.list[m], "_", direction[j], ".txt", sep=""),"\n")
#     table = read.table(file = paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/", GO.list[m], "_", direction[j], ".txt", sep=""), header=T, sep = "\t", quote="")
#     colnames(table) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
#     table.sorted = table[order(table[,5], decreasing=F),] #sorted by the degree of enrichment (GSEA) p-value
#     table.sorted = table.sorted[which(table.sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
#     table.sorted.go = table.sorted[which(table.sorted$Category == 'GOTERM_CC_FAT' |
#                                              table.sorted$Category == 'GOTERM_MF_FAT' |
#                                              table.sorted$Category == 'GOTERM_BP_FAT'),]
#     table.sorted.kegg = table.sorted[which(table.sorted$Category == 'KEGG_PATHWAY'),]
#     print(table.sorted.go[1:4,c(2,4,5,10)]) #term, percentage, pvalue, FC
#     #print(table.sorted.kegg[1:4,c(2,4,5,10)])
#     
#     #write.table(table.sorted.go, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/GO_full_", list.of.name.of.conditions[m], "_", direction[j], "_sorted.txt", sep=""), quote=F, col.names= T, row.names=F, sep = "\t")
#     #write.table(table.sorted.go[,c(2,4,5,10)], file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/GO_term_only_", list.of.name.of.conditions[m], "_", direction[j], "_sorted.txt", sep=""), quote=F, col.names= T, row.names=F, sep = "\t")
#     #write.table(table.sorted.kegg, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/KEGG_GEA_full_", list.of.name.of.conditions[m], "_", direction[j], "_sorted.txt", sep=""), quote=F, col.names= T, row.names=F, sep = "\t")
#     #write.table(table.sorted.kegg[,c(2,4,5,10)], file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_",direction[j],"regulated_DE/KEGG_GEA_term_only_", list.of.name.of.conditions[m], "_", direction[j], "_sorted.txt", sep=""), quote=F, col.names= T, row.names=F, sep = "\t")
# }
# 
