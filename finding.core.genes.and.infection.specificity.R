#Finding core genes and look into the speficifiy of infection response
#October 27th, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")
data.table <- read.table("edgeR_basic_comparison_pval_at_least_one_sig.txt",header=T)
#1. Reduce the DE dataset 
#(i) remove sterile wound and dead bacteria samples
data.table.reduced = data.table[,c(1,2,9:52)]

#(ii) Remove genes that have NON-significant p-val in all of columns
nonsig.filtered <- vector(); length.of.table <- as.numeric(dim(data.table.reduced)[1])
width.of.table <- as.numeric(dim(data.table.reduced)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...2589
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 46
        indicator <- isTRUE(data.table.reduced[k,m] > 0.05) #cutoff: FDR of 5%. True = non-significant
        if (indicator == FALSE) break  }#If the case is NOT significant, stop.
    if (m == width.of.table && indicator == T){ #if this gene has FDR>0.05 all across the conditions
        nonsig.filtered <- append(nonsig.filtered, data.table.reduced[k,1]) }
}
cat(length(nonsig.filtered), "is the number of genes that got filtered out")
filtered.data.table.reduced <- data.table.reduced[-c(nonsig.filtered),]; dim(filtered.data.table.reduced)
#Results: filtered.data.table.reduced: 2423 x 46, got rid of 166 genes.

#(iii) Convert the FDR to Y-up, Y-down, N
de.list = filtered.data.table.reduced
length.of.table <- as.numeric(dim(de.list)[1]); width.of.table <- as.numeric(dim(de.list)[2]); 
pval.indicator <- NULL; fc.indicator = NULL
for (k in 1:length.of.table){ #1, 2, 3, ... 2423
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 46
        #cat("The value for m is: ", m, "\n" )
        pval.indicator <- isTRUE(as.numeric(de.list[k,m]) > 0.05) #cutoff: FDR of 5%. True = Non-significant
        #cat("indicator: ", indicator, "\n")
        if (pval.indicator == TRUE) { #If the case is NOT significant,
            de.list[k,m] = "N"        }
        else {
            fc.indicator = isTRUE(de.list[k,m-1] >0 ) #up-regulated
            if (fc.indicator ==TRUE){
                de.list[k,m] = "Y-up"            }
            else{ de.list[k,m] = "Y-down"}
        }
    }
    pval.indicator <- NULL; fc.indicator = NULL
}

write.table(filtered.data.table.reduced, file= "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval.txt", quote=F,row.names=F)
write.table(de.list, "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", quote=F,row.names=F)

#2. Identify core genes
filtered.data.table.reduced = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval.txt", header=T)
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T)


#(i) Separate the table into 10 different conditions
M.luteus = de.list[,c(3:8)]; E.coli = de.list[,c(9:14)]; S.mar.type = de.list[,c(15:20)]; E.fae.live = de.list[,c(21:26)]; P.rett.live =de.list[,c(27:32)]
Ecc15 = de.list[,c(33:38)]; S.aureus = de.list[,c(39,40)]; P.sneebia = de.list[,c(41,42)]; S.mar.Db11 = de.list[,c(43,44)]; P.ento = de.list[,c(45,46)]
list.of.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, S.aureus, P.sneebia, S.mar.Db11, P.ento)

#(ii) Pull out genes that are up-regulated in at least one infection condition out of 10 infection conditions
sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
## Change the search word (either Y-up or Y-down) accordingly! ****************
for (i in 1:dim(de.list)[1]){ #Number of genes: 1,2,...2064
     for (j in 1:10){  #Number of conditions: 10 (M.luteus to P.entomophila)
         data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
         #print(data.to.check); print("\n")
         data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
         #print(data.to.check.for.a.gene); print("/n")
         count = length(grep("Y-down", data.to.check.for.a.gene)) ## Change the search word accordingly!
         count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
         #cat("count: ", count, "count.accumulated: ", count.accumulated); print("/n")
         #cat("j within: ",j,"  ")
         if (count > 0 ) {
             specificity=specificity + 1 #the total number of conditions that has at least Y-ups
         }
     }
    #cat("specificity: ", specificity, "\n")
    specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
    count.accumulated = 0; specificity = 0 
    if (j == 10 && count.accumulated >10){
        #print(data.to.check.for.a.gene)
        sig.in.at.least.one.conditions = rbind(sig.in.at.least.one.conditions, de.list[i,])
     }
}
library('plyr'); count(specificity.total) 
specificity.total.with.gene.name = cbind(de.list[,c(1:2)], specificity.total) #dim: 2423 x 3
int.hist = function(x,ylab="Frequency",...) {
    barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1)
}

#upregulation
specificity.total.up.reg.only = as.matrix(specificity.total.with.gene.name[which(specificity.total.with.gene.name[,3] > 0),]) #So before plotting the histogram, get rid of entries that had zero DE-up gene.
#Comments: upregulation- 1290 cases are the ones that had no significant upregulation but probably had significant downregulation.
count(specificity.total.up.reg.only[,3])
par(mfrow=c(1,1)) 
int.hist(as.numeric(specificity.total.up.reg.only[,3]), main = "Frequency of genes significantly differentially 
      upregulated in infection conditions", xlab="Number of infection conditions", ylab="Number of genes", ylim=c(0,700), xlim=c(0,10))
axis(1, at=c(1:10))

up_table = count(specificity.total.up.reg.only[,3])
barplot(up_table$freq, xlab="Number of infection conditions", ylab="Number of genes", ylim=c(0,700), xlim=c(0,12), names.arg=up_table$x, cex.axis=2.05, cex.names=2.05, cex.lab=2)


#downregulation
specificity.total.down.reg.only = as.matrix(specificity.total.with.gene.name[which(specificity.total.with.gene.name[,3] > 0),]) #So before plotting the histogram, get rid of entries that had zero DE-down gene.
#Comments: downregulation- 921 cases are the ones that had no significant upregulation but probably had significant upregulation.
count(specificity.total.down.reg.only[,3])
int.hist(as.numeric(specificity.total.down.reg.only[,3]), main = "Frequency of genes significantly differentially 
     downregulated in infection conditions", xlab="Number of infection conditions", ylab="Number of genes", ylim=c(0,700), xlim=c(0,10))
axis(1, at=c(1:10))

down_table = count(specificity.total.down.reg.only[,3])
barplot(down_table$freq, xlab="Number of infection conditions", ylab="Number of genes", ylim=c(0,700), xlim=c(0,12), names.arg=down_table$x)


#Put these two graphs together (8/2/2017)
par(mfrow = c(2,1))
par(mar=c(5,8,4,2)) #increase the margins around the graph
barplot(up_table$freq, xlab="# of infection conditions", ylab="# of upregulated genes", ylim=c(0,700), xlim=c(0,12), names.arg=up_table$x, cex.axis=1.85, cex.names=1.85, cex.lab=1.85, cex.main=2)
barplot(down_table$freq, xlab="# of infection conditions", ylab="# of downregulated genes", ylim=c(700,0), xlim=c(0,12), names.arg=down_table$x, cex.axis=1.85, cex.names=1.85, cex.lab=1.85, cex.main=2)


#(iii) List the genes with their corresponding number of infection conditions that they are positive in.
# #For upregulated:
specificity.total.up.reg.only.sorted = data.frame(specificity.total.up.reg.only[order(specificity.total.up.reg.only[,3], decreasing=T),])
colnames(specificity.total.up.reg.only.sorted) = c("gene_id", "gene_name", "number.of.sig.infection.conditions")
core.genes.10u =  specificity.total.up.reg.only.sorted[which(specificity.total.up.reg.only.sorted[,3] == "10"),]
core.genes.9u = specificity.total.up.reg.only.sorted[which(specificity.total.up.reg.only.sorted[,3] == " 9"),]
core.genes.8u = specificity.total.up.reg.only.sorted[which(specificity.total.up.reg.only.sorted[,3] == " 8"),]
write.table(specificity.total.up.reg.only.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_upregulated.txt", quote=F,row.names=F)
write.table(core.genes.10u, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_10_live_bacteria_upregulated.txt", quote=F,row.names=F)
write.table(core.genes.9u, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_9_live_bacteria_upregulated.txt", quote=F,row.names=F)
write.table(core.genes.8u, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_8_live_bacteria_upregulated.txt", quote=F,row.names=F)
write.table(count(specificity.total.up.reg.only[,3]), file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_upregulated_freq.txt", quote=F, row.names=F)

#For downregulated:
specificity.total.down.reg.only.sorted = data.frame(specificity.total.down.reg.only[order(specificity.total.down.reg.only[,3], decreasing=T),])
colnames(specificity.total.down.reg.only.sorted) = c("gene_id", "gene_name", "number.of.sig.infection.conditions")
core.genes.10d = specificity.total.down.reg.only.sorted[which(specificity.total.down.reg.only.sorted[,3] == "10"),]
core.genes.9d = specificity.total.down.reg.only.sorted[which(specificity.total.down.reg.only.sorted[,3] == " 9"),]
core.genes.8d = specificity.total.down.reg.only.sorted[which(specificity.total.down.reg.only.sorted[,3] == " 8"),]
write.table(specificity.total.down.reg.only.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_downregulated.txt", quote=F,row.names=F)
write.table(core.genes.10d, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_10_live_bacteria_downregulated.txt", quote=F,row.names=F)
write.table(core.genes.9d, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_9_live_bacteria_downregulated.txt", quote=F,row.names=F)
write.table(core.genes.8d, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_8_live_bacteria_downregulated.txt", quote=F,row.names=F)
write.table(count(specificity.total.down.reg.only[,3]), file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_downregulated_freq.txt", quote=F, row.names=F)

#(iv) Pull out gene names for GO
#UNIX: cut -f1 -d ' ' true_core_genes_all_live_bacteria_upregulated.txt > true_core_genes_live_bacteria_upregulated_gene_only.txt
#UNIX: cut -f1 -d ' ' true_core_genes_all_live_bacteria_downregulated.txt > true_core_genes_live_bacteria_downregulated_gene_only.txt

#3. Process GO analysis results -- ## Change 'up' and 'down' accordingly before starting!
#table = read.table(file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_upregulated_GO.txt", header=T, sep = "\t")
table = read.table(file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_downregulated_GO.txt", header=T, sep = "\t")

colnames(table) = c("Category","Term","Count","Percentage","PValue","Genes","List.Total","Pop.Hits","Pop.Total","Fold.Enrichment","Bonferroni","Benjamini","FDR")
table.sorted = table[order(table[,5], decreasing=F),] #sorted by the degree of enrichment (GSEA) p-value
table.sorted = table.sorted[which(table.sorted$PValue < 0.01),] #get rid of entried with EASE p-value less than 0.01
table.sorted.go = table.sorted[which(table.sorted$Category == 'GOTERM_CC_FAT' |
                                         table.sorted$Category == 'GOTERM_MF_FAT' |
                                         table.sorted$Category == 'GOTERM_BP_FAT'),]
table.sorted.kegg = table.sorted[which(table.sorted$Category == 'KEGG_PATHWAY'),]
#Upregulation
write.table(table.sorted.go, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_upregulated_GO_sorted_full.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.go[,c(2,4,5,10)], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_upregulated_GO_sorted_term_only.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.kegg, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_upregulated_KEGG_sorted_full.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.kegg[,c(2,4,5,10)], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_upregulated_KEGG_sorted_term_only.txt", quote=F, col.names= T, row.names=F, sep = "\t")
#Downregulation
write.table(table.sorted.go, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_downregulated_GO_sorted_full.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.go[,c(2,4,5,10)], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_downregulated_GO_sorted_term_only.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.kegg, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_downregulated_KEGG_sorted_full.txt", quote=F, col.names= T, row.names=F, sep = "\t")
write.table(table.sorted.kegg[,c(2,4,5,10)], file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_live_bacteria_downregulated_KEGG_sorted_term_only.txt", quote=F, col.names= T, row.names=F, sep = "\t")

###############################################################################

#4. Identify core genes by the time points #updated on 11/11/2015 from #2
rm(list=ls(all=TRUE)) #delete any previous entry
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T)

#(i) Separate the table into 10 different conditions and time points (total 22 columns)
M.luteus.12hr = de.list[,c(3:4)]; M.luteus.36hr = de.list[,c(5:6)]; M.luteus.5.5d = de.list[,c(7:8)]; E.coli.12hr = de.list[,c(9:10)]; E.coli.36hr = de.list[,c(11:12)]; E.coli.5.5d = de.list[,c(13:14)]; 
S.mar.type.12hr = de.list[,c(15:16)]; S.mar.type.36hr = de.list[,c(17:18)]; S.mar.type.5.5d = de.list[,c(19:20)]; E.fae.live.12hr = de.list[,c(21:22)]; E.fae.live.36hr = de.list[,c(23:24)]; E.fae.live.5.5d = de.list[,c(25:26)];
P.rett.live.12hr =de.list[,c(27:28)]; P.rett.live.36hr =de.list[,c(29:30)]; P.rett.live.5.5d =de.list[,c(31:32)]; Ecc15.12hr = de.list[,c(33:34)]; Ecc15.36hr = de.list[,c(35:36)]; Ecc15.5.5d = de.list[,c(37:38)]; 
S.aureus.12hr = de.list[,c(39,40)]; P.sneebia.12hr = de.list[,c(41,42)]; S.mar.Db11.12hr = de.list[,c(43,44)]; P.ento.12hr = de.list[,c(45,46)]
list.of.conditions <- list(M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d, S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d,
                           P.rett.live.12hr, P.rett.live.36hr, P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr)
list.of.name.of.conditions <- c("M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d", "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr")

#(ii) Pull out genes that are up-regulated in at least one infection condition out of 10 infection conditions
sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
spe.condition.list = vector(); spe.condition.list.total = matrix(NA, dim(de.list)[1], 1)
####################################### Change the search word (either Y-up or Y-down) accordingly!
for (i in 1:dim(de.list)[1]){ #Number of genes: 1,2,...2423
    for (j in 1:length(list.of.conditions)){  #Number of conditions: 22 (M.luteus to P.entomophila across time points)
        data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
        #print(data.to.check); print("\n")
        data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
        #print(data.to.check.for.a.gene); print("/n")
        count = length(grep("Y-down", data.to.check.for.a.gene)) ## Change the search word accordingly!
        count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
        if (count > 0 ) {
            specificity=specificity + 1 #the total number of conditions that has at least Y-ups
            spe.condition.list = paste(spe.condition.list, list.of.name.of.conditions[j], sep=", ")
        }
    }
    #cat("specificity: ", specificity, "\n")
    specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
    if (length(spe.condition.list) >0){
        spe.condition.list.total[i,1] = spe.condition.list
    }
    count.accumulated = 0; specificity = 0; spe.condition.list = vector()
#     if (j == 10 && count.accumulated >10){
#         #print(data.to.check.for.a.gene)
#         sig.in.at.least.one.conditions = rbind(sig.in.at.least.one.conditions, de.list[i,])
#     }
}

#(iii) Add the specificity number from the code #2 to the newly created listsig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
M.luteus = de.list[,c(3:8)]; E.coli = de.list[,c(9:14)]; S.mar.type = de.list[,c(15:20)]; E.fae.live = de.list[,c(21:26)]; P.rett.live =de.list[,c(27:32)]
Ecc15 = de.list[,c(33:38)]; S.aureus = de.list[,c(39,40)]; P.sneebia = de.list[,c(41,42)]; S.mar.Db11 = de.list[,c(43,44)]; P.ento = de.list[,c(45,46)]
list.of.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, S.aureus, P.sneebia, S.mar.Db11, P.ento)
sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
####################################### Change the search word (either Y-up or Y-down) accordingly!
for (i in 1:dim(de.list)[1]){ #Number of genes: 1,2,...2064
    for (j in 1:10){  #Number of conditions: 10 (M.luteus to P.entomophila)
        data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
        #print(data.to.check); print("\n")
        data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
        #print(data.to.check.for.a.gene); print("/n")
        count = length(grep("Y-down", data.to.check.for.a.gene)) ## Change the search word accordingly!
        count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
        #cat("count: ", count, "count.accumulated: ", count.accumulated); print("/n")
        #cat("j within: ",j,"  ")
        if (count > 0 ) {
            specificity=specificity + 1 #the total number of conditions that has at least Y-ups
        }
    }
    #cat("specificity: ", specificity, "\n")
    specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
    count.accumulated = 0; specificity = 0 
    if (j == 10 && count.accumulated >10){
        #print(data.to.check.for.a.gene)
        sig.in.at.least.one.conditions = rbind(sig.in.at.least.one.conditions, de.list[i,])
    }
}
specificity.total.with.gene.name = cbind(de.list[,c(1:2)], specificity.total, spe.condition.list.total) #Here, specificity.total is from the code #2
specificity.total.up.reg.only = as.matrix(specificity.total.with.gene.name[which(specificity.total.with.gene.name[,3] > 0),]) #So before plotting the histogram, get rid of entries that had zero DE-up gene (but may be DE-down)
colnames(specificity.total.up.reg.only) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "list.of.conditions.that.this.gene.was.DE")

#(iv) Add average normalized counts of selected genes
normalized.counts = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T)
normalized.counts = normalized.counts[,c(1:28)] #excluding heatkilled 
chosen.gene.ids = specificity.total.up.reg.only[,1]; count.table = data.frame()
for (i in 1:length(chosen.gene.ids)){
    index = normalized.counts[which(normalized.counts[,1] == chosen.gene.ids[i]),]
    count.table = rbind(count.table, index)
}

core.genes.with.normalized.counts = cbind(specificity.total.up.reg.only, count.table[,c(3:28)])
core.genes.with.normalized.counts.sorted = data.frame(core.genes.with.normalized.counts[order(core.genes.with.normalized.counts[,3], decreasing=T),])
write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_upregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F, row.names=F)
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_all_live_bacteria_downregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)

#5. Find genes that are induced specific to time (i.e. core time genes)
#(i) Separate the table into 10 different conditions and time points (total 22 columns)
M.luteus.12hr = de.list[,c(3:4)]; M.luteus.36hr = de.list[,c(5:6)]; M.luteus.5.5d = de.list[,c(7:8)]; E.coli.12hr = de.list[,c(9:10)]; E.coli.36hr = de.list[,c(11:12)]; E.coli.5.5d = de.list[,c(13:14)]; 
S.mar.type.12hr = de.list[,c(15:16)]; S.mar.type.36hr = de.list[,c(17:18)]; S.mar.type.5.5d = de.list[,c(19:20)]; E.fae.live.12hr = de.list[,c(21:22)]; E.fae.live.36hr = de.list[,c(23:24)]; E.fae.live.5.5d = de.list[,c(25:26)];
P.rett.live.12hr =de.list[,c(27:28)]; P.rett.live.36hr =de.list[,c(29:30)]; P.rett.live.5.5d =de.list[,c(31:32)]; Ecc15.12hr = de.list[,c(33:34)]; Ecc15.36hr = de.list[,c(35:36)]; Ecc15.5.5d = de.list[,c(37:38)]; 
S.aureus.12hr = de.list[,c(39,40)]; P.sneebia.12hr = de.list[,c(41,42)]; S.mar.Db11.12hr = de.list[,c(43,44)]; P.ento.12hr = de.list[,c(45,46)]

list.of.conditions.12hr <- list(M.luteus.12hr, E.coli.12hr,S.mar.type.12hr,E.fae.live.12hr, P.rett.live.12hr,Ecc15.12hr, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr)
list.of.conditions.36hr = list(M.luteus.36hr, E.coli.36hr,S.mar.type.36hr,E.fae.live.36hr,P.rett.live.36hr, Ecc15.36hr)
list.of.conditions.5.5d = list(M.luteus.5.5d, E.coli.5.5d, S.mar.type.5.5d, E.fae.live.5.5d, P.rett.live.5.5d, Ecc15.5.5d)
list.of.name.of.conditions.12hr <- c("M.luteus.12hr","E.coli.12hr","S.mar.type.12hr", "E.fae.live.12hr","P.rett.live.12hr", "Ecc15.12hr","S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr")
list.of.name.of.conditions.36hr <- c("M.luteus.36hr", "E.coli.36hr","S.mar.type.36hr", "E.fae.live.36hr","P.rett.live.36hr", "Ecc15.36hr")
list.of.name.of.conditions.5.5d <- c("M.luteus.5.5d",   "E.coli.5.5d",  "S.mar.type.5.5d",  "E.fae.live.5.5d", "P.rett.live.5.5d",  "Ecc15.5.5d")

#(ii) Pull out genes that are up-regulated in at least one infection condition out of 10 infection conditions
sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
spe.condition.list = vector(); spe.condition.list.total = matrix(NA, dim(de.list)[1], 1)
####################################### 
#1. Change the search word (either Y-up or Y-down) accordingly!
#2. Plug in the right list of conditions and the list of name of conditions
list.of.conditions = list.of.conditions.5.5d
list.of.name.of.conditions =list.of.name.of.conditions.5.5d
for (i in 1:dim(de.list)[1]){ #Number of genes: 1,2,...2423
    for (j in 1:length(list.of.conditions)){  #Number of conditions: 10 if 12hr, 6 if 36hr
        data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
        data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
        count = length(grep("Y-down", data.to.check.for.a.gene)) ## Change the search word accordingly!
        count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
        if (count > 0 ) {
            specificity=specificity + 1 #the total number of conditions that has at least Y-ups
            spe.condition.list = paste(spe.condition.list, list.of.name.of.conditions[j], sep=", ")
        }
    }
    #cat("specificity: ", specificity, "\n")
    specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
    if (length(spe.condition.list) >0){
        spe.condition.list.total[i,1] = spe.condition.list
    }
    count.accumulated = 0; specificity = 0; spe.condition.list = vector()
}

#(iii) Put together the gene id, name, number of specificity cases, and the list of conditions
specificity.total.with.gene.name = cbind(de.list[,c(1:2)], specificity.total, spe.condition.list.total) 
specificity.total.reg.only = as.matrix(specificity.total.with.gene.name[which(specificity.total.with.gene.name[,3] > 0),]) #So before plotting the histogram, get rid of entries that had zero DE-up gene (but may be DE-down)
colnames(specificity.total.reg.only) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "list.of.conditions.that.this.gene.was.DE")

#(iv) Add average normalized counts of selected genes
normalized.counts = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T)
normalized.counts = normalized.counts[,c(1:28)] #excluding heatkilled 
chosen.gene.ids = specificity.total.reg.only[,1]; count.table = data.frame()
for (i in 1:length(chosen.gene.ids)){
    index = normalized.counts[which(normalized.counts[,1] == chosen.gene.ids[i]),]
    count.table = rbind(count.table, index)
}

core.genes.with.normalized.counts = cbind(specificity.total.reg.only, count.table[,c(3:28)])
core.genes.with.normalized.counts.sorted = data.frame(core.genes.with.normalized.counts[order(core.genes.with.normalized.counts[,3], decreasing=T),])
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_12hr_all_live_bacteria_upregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_12hr_all_live_bacteria_downregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_36hr_all_live_bacteria_upregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_36hr_all_live_bacteria_downregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)
#write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_5.5d_all_live_bacteria_upregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)
write.table(core.genes.with.normalized.counts.sorted, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_time_points/core_genes_by_5.5d_all_live_bacteria_downregulated_with_conditions_and_normalized_counts.txt", sep="\t", quote=F,row.names=F)

#6. Compare the "core genes" to sterile.wound DE genes in UNIX = #ontrast the core up/down-regulated genes to sterile.wound DE genes
#If the core gene shows up as DEG in any time point of sterile wound, it was excluded
###############################################################################

#i) Comparing upDEG-core.8-10 to the accumulated list of sterile wound across time points (clean.prick)
#cd /Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs
#cat list_of_upregulated_DE_genes_in_clean.prick.12hr.txt list_of_upregulated_DE_genes_in_clean.prick.36hr.txt list_of_upregulated_DE_genes_in_clean.prick.5.5d.txt > list_of_upregulated_DE_genes_in_clean.prick.txt
#cut -f2 -d " "  list_of_upregulated_DE_genes_in_clean.prick.txt > list_of_upregulated_DE_genes_in_clean.prick_name_only.txt
clean.prick = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs/list_of_upregulated_DE_genes_in_clean.prick_name_only.txt", header=T)
clean.prick.uniq = as.matrix(unique(clean.prick)) #87
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_all_live_bacteria_upregulated.txt", header=T)
core.genes.names = matrix(core.genes[,2]) #135 genes
non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = clean.prick.uniq[which(clean.prick.uniq[,1] == core.genes.names[i,1]),]
    if (length(index)<1){ #if there is no overlap
        non.overlap = rbind(non.overlap,core.genes.names[i,1])
    }
}
dim(non.overlap) #shows the number of reduced core genes that do not have overlap with the sterile wound = 71 genes
write.table(non.overlap, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_live_bacteria_upregulated_clean_prick_DEGs_removed.txt", quote=F, row.names=F)

#ii) Comparing downDEG-core.8-10 to the accumulated list of sterile wound across time points (clean.prick)
#cd /Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs_in_each_conditions
#cat list_of_downregulated_DE_genes_in_clean.prick.12hr.txt list_of_downregulated_DE_genes_in_clean.prick.36hr.txt list_of_downregulated_DE_genes_in_clean.prick.5.5d.txt > list_of_downregulated_DE_genes_in_clean.prick.txt
#cut -f2 -d " "  list_of_downregulated_DE_genes_in_clean.prick.txt > list_of_downregulated_DE_genes_in_clean.prick_name_only.txt
clean.prick = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/list_of_DEGs/list_of_downregulated_DE_genes_in_clean.prick_name_only.txt", header=T)
clean.prick.uniq = as.matrix(unique(clean.prick)) #29 genes
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_all_live_bacteria_downregulated.txt", header=T)
core.genes.names = matrix(core.genes[,2]) #54 genes
non.overlap = vector()
for (i in 1:dim(core.genes.names)[1]){
    index = clean.prick.uniq[which(clean.prick.uniq[,1] == core.genes.names[i,1]),]
    if (length(index)<1){ #if there is no overlap
        non.overlap = rbind(non.overlap,core.genes.names[i,1])
    }
}
dim(non.overlap) #shows the number of reduced core genes that do not have overlap with the sterile wound = 42 genes
write.table(non.overlap, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_live_bacteria_downregulated_clean_prick_DEGs_removed.txt", quote=F, row.names=F)


#7. Compare the list of core genes to the DIRGs-385 genes (upregulation only)
###############################################################################
#i) Compare the list of true core genes (8+) to the 385 DIRGs
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_all_live_bacteria_upregulated.txt", header=T)
core.genes.names = matrix(core.genes[,1]) #135 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_385_updated_simple.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
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
dim(overlap) #shows the genes that overlap between my true core genes and DIRGs = 74 genes
dim(non.overlap) #shows the genes that overlap between my true core genes and DIRGs = 61 genes
overlap.2 = cbind(overlap,c(rep("Y",dim(overlap)[1]))); non.overlap.2 = cbind(non.overlap,c(rep("N",dim(non.overlap)[1])))
colnames(overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "DIRGs")
colnames(non.overlap.2) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "DIRGs")
core.genes.with.dirg = rbind(overlap.2, non.overlap.2); core.genes.with.dirg.sorted = core.genes.with.dirg[order(as.numeric(rownames(core.genes.with.dirg)), decreasing=F),]
write.table(core.genes.with.dirg.sorted, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/true_core_genes_all_live_bacteria_upregulated_compared_to_DIRGs.txt", quote=F, row.names=F, col.names=T)

#ii) Compare the list of all DEGs (1+ in at least one infection, including clean prick, infected, and heatkilled) to the 385 DIRGs
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
core.genes.names = matrix(core.genes[,1]) #2589 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_385_updated_simple.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
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
dim(overlap) #shows the genes that overlap between my true core genes and DIRGs = 228 genes
dim(non.overlap) #shows the genes that overlap between my true core genes and DIRGs = 2361 genes
yes = matrix(c(rep("Y",dim(overlap)[1]))); no = matrix(c(rep("N",dim(non.overlap)[1]))); colnames(yes) = c("DIRGs"); colnames(no) = c("DIRGs")
overlap.2 = cbind(overlap[,c(1:2)], yes, overlap[,c(3:64)]); non.overlap.2 = cbind(non.overlap[,c(1:2)], no, non.overlap[,c(3:64)])
core.genes.with.dirg = rbind(overlap.2, non.overlap.2)
core.genes.with.dirg.sorted = core.genes.with.dirg[order(as.numeric(rownames(core.genes.with.dirg)), decreasing=F),]
write.table(core.genes.with.dirg.sorted, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/DEG-ups_in_all_conditions_compared_to_DIRGs.txt", quote=F, row.names=F, col.names=T)

#iii) Compare the list of all DEGs (1+ in at least one live infection, excluding clean prick, infected, and heatkilled) to the 385 DIRGs
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
core.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/core_genes_by_num_of_conditions/core_genes_all_live_bacteria_upregulated.txt", header=T)
core.genes.names = matrix(core.genes[,1]) #1286 genes
dirg = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(8_or_more_conditions)/DIRG_385_updated_simple.txt", header=T) #The original CG numbers are converted into FBgn Flybase IDs (11/12/2015)
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


#8. Determine if any of the TFs, kinase-phosph, serine proteases are considered significant in any conditions
###############################################################################
#i) Read in the genes and make comparisons
rm(list=ls(all=TRUE)) 
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/")
live.bacteria.table = read.table("basic_comparison_reduced_to_live_bacteria_samples_pval_let_converted.txt", header=T)
tfs = read.table("tfs_to_test.txt", header =T) #715 genes
kinase = read.table("kinase_phosph_to_test.txt", header =T) #288 genes
serine = read.csv("clip-domain_serine_proteases_to_test_Veillard_et_al.txt", header=T) #42 genes
functional.genes.type = c("TFs", "Kinase", "Serine.proteases"); functional.genes = list(tfs, kinase, serine)
#ii) Organize the list of DE-TFs, DE-Kinase, or DE-Serine proteases in the number and the type of conditions in which they were considered as DE
for (k in 1:length(functional.genes)){
    data.to.work.with = data.frame(functional.genes[k]); overlap = vector()
    for (i in 1:dim(live.bacteria.table)[1]){
        index = live.bacteria.table[which(live.bacteria.table[,1] == as.character(data.to.work.with[i,1])),]; #print(index)
        if (dim(index)[1] >= 1){ #if they overlap
            overlap = rbind(overlap, index)
        }
    }
    print(dim(overlap)) #70 TF are DE in some conditions; 40 Kinases are DE in some conditions; 16 serine-protease are DE in some conditions
    write.table(overlap, file=paste(functional.genes.type[k],"_in_live_bacteria_samples_pval_let_converted.txt",sep=""), quote=F, row.names=F, col.names=T)
    
    M.luteus.12hr = overlap[,c(3:4)]; M.luteus.36hr = overlap[,c(5:6)]; M.luteus.5.5d = overlap[,c(7:8)]; E.coli.12hr = overlap[,c(9:10)]; E.coli.36hr = overlap[,c(11:12)]; E.coli.5.5d = overlap[,c(13:14)]; 
    S.mar.type.12hr = overlap[,c(15:16)]; S.mar.type.36hr = overlap[,c(17:18)]; S.mar.type.5.5d = overlap[,c(19:20)]; E.fae.live.12hr = overlap[,c(21:22)]; E.fae.live.36hr = overlap[,c(23:24)]; E.fae.live.5.5d = overlap[,c(25:26)];
    P.rett.live.12hr =overlap[,c(27:28)]; P.rett.live.36hr =overlap[,c(29:30)]; P.rett.live.5.5d =overlap[,c(31:32)]; Ecc15.12hr = overlap[,c(33:34)]; Ecc15.36hr = overlap[,c(35:36)]; Ecc15.5.5d = overlap[,c(37:38)]; 
    S.aureus.12hr = overlap[,c(39,40)]; P.sneebia.12hr = overlap[,c(41,42)]; S.mar.Db11.12hr = overlap[,c(43,44)]; P.ento.12hr = overlap[,c(45,46)]
    list.of.conditions <- list(M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d, S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d,
                               P.rett.live.12hr, P.rett.live.36hr, P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr)
    list.of.name.of.conditions <- c("M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d", "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                    "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr")
    
    sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector() #Pull out genes that are up-regulated in at least one infection condition out of 10 infection conditions
    spe.condition.list = vector(); spe.condition.list.total = matrix(NA, dim(overlap)[1], 1)
    ####################################### Change the search word (either Y-up or Y-down) accordingly!
    for (i in 1:dim(overlap)[1]){ #Number of genes: 1,2,...n
        for (j in 1:length(list.of.conditions)){  #Number of conditions: 22 (M.luteus to P.entomophila across time points)
            data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
            #print(data.to.check); print("\n")
            data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
            #print(data.to.check.for.a.gene); print("/n")
            count = length(grep("Y-up", data.to.check.for.a.gene)) ## Change the search word accordingly!
            count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
            if (count > 0 ) {
                specificity=specificity + 1 #the total number of conditions that has at least Y-ups
                spe.condition.list = paste(spe.condition.list, list.of.name.of.conditions[j], sep=", ")
            }
        }
        #cat("specificity: ", specificity, "\n")
        specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
        if (length(spe.condition.list) >0){
            spe.condition.list.total[i,1] = spe.condition.list
        }
        count.accumulated = 0; specificity = 0; spe.condition.list = vector()
        #     if (j == 10 && count.accumulated >10){
        #         #print(data.to.check.for.a.gene)
        #         sig.in.at.least.one.conditions = rbind(sig.in.at.least.one.conditions, de.list[i,])
        #     }
    }
    
    #Add the specificity number from the code #2 to the newly created listsig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
    M.luteus = overlap[,c(3:8)]; E.coli = overlap[,c(9:14)]; S.mar.type = overlap[,c(15:20)]; E.fae.live = overlap[,c(21:26)]; P.rett.live =overlap[,c(27:32)]
    Ecc15 = overlap[,c(33:38)]; S.aureus = overlap[,c(39,40)]; P.sneebia = overlap[,c(41,42)]; S.mar.Db11 = overlap[,c(43,44)]; P.ento = overlap[,c(45,46)]
    list.of.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, S.aureus, P.sneebia, S.mar.Db11, P.ento)
    sig.in.at.least.one.conditions = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
    ####################################### Change the search word (either Y-up or Y-down) accordingly!
    for (i in 1:dim(overlap)[1]){ #Number of genes: 1,2,...n
        for (j in 1:10){  #Number of conditions: 10 (M.luteus to P.entomophila)
            data.to.check = as.matrix(as.data.frame(list.of.conditions[j]))
            #print(data.to.check); print("\n")
            data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
            #print(data.to.check.for.a.gene); print("/n")
            count = length(grep("Y-up", data.to.check.for.a.gene)) ## Change the search word accordingly!
            count.accumulated = count.accumulated + count #the total number of Y-ups in a given gene (regardless of time points)
            #cat("count: ", count, "count.accumulated: ", count.accumulated); print("/n")
            #cat("j within: ",j,"  ")
            if (count > 0 ) {
                specificity=specificity + 1 #the total number of conditions that has at least Y-ups
            }
        }
        #cat("specificity: ", specificity, "\n")
        specificity.total = rbind(specificity.total,specificity) #A vector for a histogram!!!!
        count.accumulated = 0; specificity = 0 
        if (j == 10 && count.accumulated >10){
            #print(data.to.check.for.a.gene)
            sig.in.at.least.one.conditions = rbind(sig.in.at.least.one.conditions, de.list[i,])
        }
    }
    specificity.total.with.gene.name = cbind(overlap[,c(1:2)], specificity.total, spe.condition.list.total) #Here, specificity.total is from the code #2
    specificity.total.reg.only = as.matrix(specificity.total.with.gene.name[which(specificity.total.with.gene.name[,3] > 0),]) #So before plotting the histogram, get rid of entries that had zero DE-up gene (but may be DE-down)
    colnames(specificity.total.reg.only) = c("gene_id", "gene_name", "number.of.sig.infection.conditions", "list.of.conditions.that.this.gene.was.DE")
    specificity.total.reg.only.sorted = data.frame(specificity.total.reg.only[order(specificity.total.reg.only[,3], decreasing=T),])
    write.table(specificity.total.reg.only.sorted, file = paste(functional.genes.type[k],"_DE-ups_in_live_bacteria_samples.txt",sep=""), quote=F, row.names=F )
    
}
