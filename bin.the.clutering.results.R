#Bin the clustering results
#August 19, 2016 and updated on August 23, 2016
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/EBSeq-HMM_Aug_2016/")

#1. Create a list of genes from 8 conditions and only pick the unique entries
M.luteus = read.table("M.luteus.with.most.liekly.paths.based.on.normalized.counts.txt", header=T); E.coli = read.table("E.coli.with.most.liekly.paths.based.on.normalized.counts.txt", header=T)
S.mar.type = read.table("S.mar.type.with.most.liekly.paths.based.on.normalized.counts.txt", header=T); Ecc15 = read.table("Ecc15.with.most.liekly.paths.based.on.normalized.counts.txt", header=T)
P.rett.live = read.table("P.rett.live.with.most.liekly.paths.based.on.normalized.counts.txt", header=T); E.fae.live = read.table("E.fae.live.with.most.liekly.paths.based.on.normalized.counts.txt", header=T)
P.rett.heatkilled = read.table("P.rett.heatkilled.with.most.liekly.paths.based.on.normalized.counts.txt", header=T); E.fae.heatkilled = read.table("E.fae.heatkilled.with.most.liekly.paths.based.on.normalized.counts.txt", header=T)
gene.list = c(row.names(M.luteus), row.names(E.coli), row.names(S.mar.type), row.names(Ecc15), row.names(P.rett.live), row.names(E.fae.live), row.names(P.rett.heatkilled), row.names(E.fae.heatkilled))
gene.list.uniq = unique(gene.list) #951 genes (confident calls only)

#Put gene name along with FBgn number
full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
output.mx <- matrix(NA, ncol=1, nrow=length(gene.list.uniq)); colnames(output.mx) <-c("gene_name")
for (i in 1:length(gene.list.uniq)){
    gene.to.compare = as.character(gene.list.uniq[i])
    checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
    output.mx[i,1] <- as.character(checking.against.full.list[1,10])
}

#2. Create a table that shows the clustering pattern information
cluster.table = matrix(NA, nrow=length(gene.list.uniq), ncol=10); cluster.table[,1] = gene.list.uniq; cluster.table[,10] = output.mx
colnames(cluster.table) = c("gene.id", "M.luteus","E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled","gene.name")
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

for (i in 1:length(list.of.bacteria.name)){
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[i])); count.table <- condition.count.table
    for (j in 1:dim(count.table)[1]){ #for a given condition
        gene.name = row.names(count.table)[j]
        pattern = count.table[j,1]
        cluster.table[which(cluster.table[,1] == gene.name),][i+1] = pattern
    }
}
write.table(cluster.table, file="analysis.Aug23/bin.table.of.the.cluster.results.across.infection.conditions.txt", quote=F, row.names = F, col.names = T)


#3. What are the genes that behave differently when infecting bacteria is live or heatkilled?
cluster.table.E.fae = cluster.table[,c(1,5,8,10)]; cluster.table.E.fae.na.omit = na.omit(cluster.table.E.fae) #74 genes
cluster.table.P.rett = cluster.table[,c(1,6,9,10)]; cluster.table.P.rett.na.omit = na.omit(cluster.table.P.rett) #33 gens

#E.faecalis live vs E.faecalis heatkilled. Out of 74 genes, 7 genes have different patterns between conditions
E.faecalis.live.vs.heatkilled.infection.diff.genes = c()
for (k in 1:dim(cluster.table.E.fae.na.omit)[1]){
    if (cluster.table.E.fae.na.omit[k,2] != cluster.table.E.fae.na.omit[k,3]) { #If the gene patterns are not the same
        E.faecalis.live.vs.heatkilled.infection.diff.genes = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, cluster.table.E.fae.na.omit[k,])
    }
}
write.table(E.faecalis.live.vs.heatkilled.infection.diff.genes, file="analysis.Aug23/genes.that.behave.differently.between.E.faecalis.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#P.rettgeri live vs P.rettgeri heatkilled. Out of 33 genes, 11 genes have different patterns between conditions
P.rettgeri.live.vs.heatkilled.infection.diff.genes = c()
for (n in 1:dim(cluster.table.P.rett.na.omit)[1]){
    if (cluster.table.P.rett.na.omit[n,2] != cluster.table.P.rett.na.omit[n,3]) { #If the gene patterns are not the same
        P.rettgeri.live.vs.heatkilled.infection.diff.genes = rbind(P.rettgeri.live.vs.heatkilled.infection.diff.genes, cluster.table.P.rett.na.omit[n,])
    }
}
write.table(P.rettgeri.live.vs.heatkilled.infection.diff.genes, file="analysis.Aug23/genes.that.behave.differently.between.P.rettgeri.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#Are there common patterns of live vs heatkilled bacteria? #No, there was no gene that had a common pattern between two conditions (this doesn't mean that there is no gene. It's simply missing in data)
cluster.table.heat = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, P.rettgeri.live.vs.heatkilled.infection.diff.genes)
dup_idx = duplicated(cluster.table.heat[,1]); cluster.table.heat.with.idx = cbind(cluster.table.heat, dup_idx)
cluster.table.heat.common = cluster.table.heat.with.idx[which(cluster.table.heat.with.idx[,5] == "TRUE"),]; cluster.table.heat.common = c(cluster.table.heat.common[,1])

cluster.table.heat.common.table = c()
for(a in 1:length(cluster.table.heat.common)){
    line = cluster.table[which(cluster.table[,1] == cluster.table.heat.common[a]),]
    cluster.table.heat.common.table = rbind(cluster.table.heat.common.table, line)
}
write.table(cluster.table.heat.common.table[,c(1,5,8,6,9,10)], file = "analysis.Aug23/genes.that.behave.differently.between.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)
#Answer: None.

#4. What are the genes that behave the same way across infection conditions? ("core of core genes")
cluster.table.live = cluster.table[,c(1:7,10)]; cluster.table.live.na.omit = na.omit(cluster.table.live) #19 genes that have information across conditions
core.of.core=c()
for (b in 1:dim(cluster.table.live.na.omit)[1]){
    first.pattern = as.character(cluster.table.live.na.omit [b,2])
    if (first.pattern == cluster.table.live.na.omit [b,3] && first.pattern == cluster.table.live.na.omit [b,4] && first.pattern == cluster.table.live.na.omit [b,5] && first.pattern == cluster.table.live.na.omit [b,6] && first.pattern == cluster.table.live.na.omit [b,7]){
        core.of.core = rbind(core.of.core, cluster.table.live.na.omit [b,])
    }
}
write.table(core.of.core, file = "analysis.Aug23/genes.that.behave.the.same.across.infection.conditions.txt", quote=F, col.names=T, row.names = F) #16 genes


#5. What are the genes that behave the same way in the Gram-positive infection (M.luteus, E.fae.live)?
cluster.table.positive = cluster.table[,c(1,2,5,10)]; cluster.table.positive.na.omit = na.omit(cluster.table.positive) #155 genes that have information across two Gram positive infections
gram.positive.only=c()
for (c in 1:dim(cluster.table.positive.na.omit)[1]){
    first.pattern = as.character(cluster.table.positive.na.omit[c,2]) #M.luteus pattern
    if (first.pattern == cluster.table.positive.na.omit[c,3]){
        gram.positive.only = rbind(gram.positive.only, cluster.table.positive.na.omit[c,])
    }
}
write.table(gram.positive.only, file = "analysis.Aug23/genes.that.behave.the.same.in.gram.positive.infections.txt", quote=F, col.names=T, row.names = F) #150 genes


#6. What are the genes that behave the same way in the Gram-negative infection (E.coli, S.marcescens Type, P.rettgeri, Ecc15)?
cluster.table.negative = cluster.table[,c(1,3,4,6,7,10)]; cluster.table.negative.na.omit = na.omit(cluster.table.negative) #14 genes that have information across two Gram positive infections
gram.negative.only=c()
for (d in 1:dim(cluster.table.negative.na.omit)[1]){
    first.pattern = as.character(cluster.table.negative.na.omit[d,2]) #E.coli pattern
    if (first.pattern == cluster.table.negative.na.omit[d,3] && first.pattern == cluster.table.negative.na.omit[d,4] && first.pattern == cluster.table.negative.na.omit[d,5]){
        gram.negative.only = rbind(gram.negative.only, cluster.table.negative.na.omit[d,])
    }
}
#All the core genes are a part of the Gram negative only genes. In other words, these genes behave the same way all across the infection conditions.
write.table(gram.negative.only, file = "analysis.Aug23/genes.that.behave.the.same.in.gram.negative.infections.txt", quote=F, col.names=T, row.names = F) #18 genes


#7. Check the expression of individual genes
rnaseq.table = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T)
rnaseq.table[which(rnaseq.table$gene.id == "FBgn0002563"),][,3]

#Is there a gene that behaves differently depending on infection conditions?
cluster.table.live.na.omit #23 genes. No, there's no gene that has a different pattern in each infection condition.
