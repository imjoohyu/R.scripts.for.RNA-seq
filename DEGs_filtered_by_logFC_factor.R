#DEGs filtered by logFC factor
#March 17th, 2017
#Joo Hyun Im (ji72)


rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")

DEGs_total = read.table("edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 x 64


#1 Does the number of DEGs change when the more strict filter is applied?
list.of.name.of.conditions <- c("Clean.prick.12hr", "Clean.prick.36hr", "Clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d","S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr","P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr","E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
count=1

table=data.frame()
for (n in seq(4,64,2)){
    sig.genes.in.this.condition = subset(DEGs_total, grepl("Y",DEGs_total[,n]) )
    subset = sig.genes.in.this.condition[,c(1:2,n-1,n)]
    subset.sorted = subset[order(subset[,3], decreasing=T),]

        #write.table(subset.sorted[which(subset.sorted[,3] > 1),], file=paste("../list_of_DEGs_factor_of_two_or_more/list_of_upregulated_DE_genes_FC1_more_in_",list.of.name.of.conditions[count],".txt",sep=""), quote=F, col.names= T, row.names=F)
    #write.table(subset.sorted[which(subset.sorted[,3] < -1),], file=paste("../list_of_DEGs_factor_of_two_or_more/list_of_downregulated_DE_genes_FC1_more_in_",list.of.name.of.conditions[count],".txt",sep=""), quote=F, col.names= T, row.names=F)
    #cat("# of up-DEGs for", list.of.name.of.conditions[count], " is ", dim(subset.sorted[which(subset.sorted[,3] > 1),])[1], "\n")
    #cat("# of down-DEGs for", list.of.name.of.conditions[count], " is ", dim(subset.sorted[which(subset.sorted[,3] < -1),])[1], "\n")
    
    table_row = cbind(list.of.name.of.conditions[count], dim(subset.sorted[which(subset.sorted[,3] > 0),])[1], dim(subset.sorted[which(subset.sorted[,3] < 0),])[1], dim(subset.sorted[which(subset.sorted[,3] > 1),])[1], dim(subset.sorted[which(subset.sorted[,3] < -1),])[1])
    table = rbind(table, table_row)
    count = count +1
}
colnames(table) = c("condition", "num.of.up.DEGs", "num.of.down.DEGs", "num.of.up.DEGs.FC1", "num.of.down.DEGs.FC1")
write.table(table, file="../list_of_DEGs_factor_of_two_or_more/comparison_of_DEG_number_with_and_without_FC1_filter.txt", quote=F, row.names = F, col.names = T)

#2 Does the composition of the core genes drastically change when the more strict filter is applied?




