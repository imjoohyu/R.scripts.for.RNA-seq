#Findinh chronic-specific genes
#December 7, 2015
#Joo Hyun Im (ji72)

#Find chronic-infection specific genes: Are there genes that are not DE in 12hr and 36hr, but are DE in 5.5d?

rm(list=ls(all=TRUE)) #delete any previous entry
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T) #2423 genes

#1. Choose genes that have non-significant p-values in 12hr and 36hr
filter <- vector(); indicator <- NULL
for (i in 1:dim(de.list)[1]){ # each gene from 2423 genes (from live infections)
    for (j in c(4,6,10,12,16,18,22,24,28,30,34,36,40,42,44,46)){ #targetting 12hr and 36hr
        indicator <- isTRUE(de.list[i,j] == "N" ) #Indicate TRUE if this gene i is not DE in condition j
        #cat("j: ", j," indicator: ",indicator)
        if (indicator == FALSE) break  }#If the gene is significant, stop.
    if (j == 46 && indicator == TRUE){ #if this gene i is not DE in all conditions
        filter <- append(filter, de.list[i,1])
    }
}
cat(length(filter), "is the number of genes that are not DE in 12hr AND 36hr") #238 genes
de.list.reduced <- de.list[c(filter),]; dim(de.list.reduced)

# #2. Choose genes that have significant p-values in 5.5d -- no need to do this because the results from #1 are up or downregulated in at least one 5.5d conditions 
# filter <- vector(); indicator <- NULL
# for (i in 1:dim(de.list.reduced)[1]){ # each gene from 238 genes (from live infections)
#     for (j in c(8,14,20,26,32,38)){ #targetting 5.5d for M.leus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15
#         indicator <- isTRUE(de.list.reduced[i,j] == "N" ) #Indicate TRUE if this gene i is not DE in condition j
#         if (indicator == FALSE) break  }#If the gene is significant, stop.
#     if (j == 46 && indicator == TRUE){ #if this gene i is not DE in all conditions
#         filter <- append(filter, de.list.reduced[i,1])
#     }
# }
# cat(length(filter), "is the number of genes that are not DE in 12hr AND 36hr BUT are DE in any of the 5.5d samples") #0 genes 
write.table(de.list.reduced, "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/chronic-infection-specific.genes_all_time_points.txt", quote=F, row.names=F)
write.table(de.list.reduced[,c(1,2,7,8,13,14,19,20,25,26,31,32,37,38)], "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/chronic-infection-specific.genes_5.5d_only.txt", quote=F, row.names=F)

#3. Split the 5.5d-only-DE genes by upregulated DEGs and downregulated DEGs
filter <- vector(); filter.up = vector(); filter.down = vector(); indicator <- NULL; direction = c("Y-up", "Y-down")
rownames(de.list.reduced) = c(1:238); head(de.list.reduced)
for (p in 1:2){
    for (m in 1:dim(de.list.reduced)[1]){
        for (j in c(8,14,20,26,32,38)){ #targetting 5.5d for M.leus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15
            #print(de.list.reduced[m,j])
            indicator <- isTRUE(de.list.reduced[m,j] == direction[p]) #Indicate TRUE if this gene i is [direction]-DEG in condition j
            #cat("m: ", i, "j: ", j, " indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the gene is significant, add it to the filter, and then stop.
                filter <- append(filter, rownames(de.list.reduced[m,]))}
            if (indicator == TRUE) break}
    }
    if(p == 1) {filter.up = filter} #73 genes
    if(p == 2) {filter.down = filter} #165 genes
    filter = vector()
}
de.list.reduced.deg.up <- de.list.reduced[c(filter.up),]; de.list.reduced.deg.down <- de.list.reduced[c(filter.down),]
write.table(de.list.reduced.deg.up[,c(1,2,7,8,13,14,19,20,25,26,31,32,37,38)], "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/chronic-infection-specific.upregulated_genes_5.5d_only.txt", quote=F, row.names=F)
write.table(de.list.reduced.deg.down[,c(1,2,7,8,13,14,19,20,25,26,31,32,37,38)], "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/chronic-infection-specific.downregulated_genes_5.5d_only.txt", quote=F, row.names=F)

#Comments (12/08/2015): 
#GO term analysis results indicate that the upregulated genes' common (3+) GO term is:
#'transferase' ("Enzyme that transfers a chemical group, e.g. a methyl group or a glycosyl group from one compound (donor) to another compound (acceptor).")
#GO term analysis results indicate that the downregulated genes' common (3+) GO term is:
#None

