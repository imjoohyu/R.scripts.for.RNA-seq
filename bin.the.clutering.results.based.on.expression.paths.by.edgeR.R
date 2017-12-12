#Bin the clustering results based on edgeR expression path assignment
#October 12, 2016
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

########################################################################################################
#Part I. Unchallenged vs Infected dataset
########################################################################################################

#1. What are the genes that behave differently when infecting bacteria is live or heatkilled?
cluster.table = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_NS_removed.txt", header=T)
cluster.table.E.fae = cluster.table[,c(1,2,6,9)]; cluster.table.P.rett = cluster.table[,c(1,2,7,10)]
cluster.table.E.fae$E.fae.live.path = factor(cluster.table.E.fae$E.fae.live.path); cluster.table.E.fae$E.fae.heatkilled.path = factor(cluster.table.E.fae$E.fae.heatkilled.path)
cluster.table.P.rett$P.rett.live.path = factor(cluster.table.P.rett$P.rett.live.path); cluster.table.P.rett$P.rett.heatkilled.path = factor(cluster.table.P.rett$P.rett.heatkilled.path)

#E.faecalis live vs E.faecalis heatkilled. 884 genes have different patterns between conditions
E.faecalis.live.vs.heatkilled.infection.diff.genes = c()
for (k in 1:dim(cluster.table.E.fae)[1]){
    if (as.character(cluster.table.E.fae[k,3]) != as.character(cluster.table.E.fae[k,4])) { #If the gene patterns are not the same
        E.faecalis.live.vs.heatkilled.infection.diff.genes = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, cluster.table.E.fae[k,])
    }
}
write.table(E.faecalis.live.vs.heatkilled.infection.diff.genes, file="binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.differently.between.E.faecalis.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#P.rettgeri live vs P.rettgeri heatkilled. 1334 genes have different patterns between conditions
P.rettgeri.live.vs.heatkilled.infection.diff.genes = c()
for (n in 1:dim(cluster.table.P.rett)[1]){
    if (as.character(cluster.table.P.rett[n,3]) != as.character(cluster.table.P.rett[n,4])) { #If the gene patterns are not the same
        P.rettgeri.live.vs.heatkilled.infection.diff.genes = rbind(P.rettgeri.live.vs.heatkilled.infection.diff.genes, cluster.table.P.rett[n,])
    }
}
write.table(P.rettgeri.live.vs.heatkilled.infection.diff.genes, file="binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.differently.between.P.rettgeri.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#2. Are there common patterns of live vs heatkilled bacteria?
#Part 1. Pick genes that have changed between live and heatkilled in both infections.
colnames(E.faecalis.live.vs.heatkilled.infection.diff.genes)[3:4]= c("live","heatkilled")
colnames(P.rettgeri.live.vs.heatkilled.infection.diff.genes)[3:4]= c("live","heatkilled") #unify the column names
cluster.table.heat = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, P.rettgeri.live.vs.heatkilled.infection.diff.genes)
dup_idx = duplicated(cluster.table.heat[,1]); cluster.table.heat.with.idx = cbind(cluster.table.heat, dup_idx) #pick genes that were different between live and heatkilled conditions in each infection and that had information in both infections
cluster.table.heat.common = cluster.table.heat.with.idx[which(cluster.table.heat.with.idx$dup_idx == "TRUE"),]; cluster.table.heat.common = c(as.character(cluster.table.heat.common[,1]))

cluster.table.heat.common.table = c()
for(a in 1:length(cluster.table.heat.common)){
    line = cluster.table[which(cluster.table[,1] == cluster.table.heat.common[a]),]
    cluster.table.heat.common.table = rbind(cluster.table.heat.common.table, line)
}
#Part 2. Find out which genes have the same pattern between live and heatkilled in both bacterial conditions.
cluster.table.heat.common.table = cluster.table.heat.common.table[,c(1,2,6,9,7,10)]
cluster.table.heat.same.table = c()
for (b in 1:dim(cluster.table.heat.common.table)[1]){
    if (as.character(cluster.table.heat.common.table[b,3]) == as.character(cluster.table.heat.common.table[b,5])
        && as.character(cluster.table.heat.common.table[b,4]) == as.character(cluster.table.heat.common.table[b,6])
        ) { #If the gene patterns are not the same
        cluster.table.heat.same.table = rbind(cluster.table.heat.same.table, cluster.table.heat.common.table[b,])
    }
}

write.table(cluster.table.heat.same.table, file = "binning.the.clustering.results/edgeR_unchallenged_vs_infected_common_genes.that.behave.differently.between.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)
#Answer: 194 genes.


#3. What are the genes that behave the same way across infection conditions? ("core of core genes")
#Part 1. Pick out genes that have the same expression path across all conditions
cluster.table.live = cluster.table[,c(1:8)]; core.of.core=c()
for (b in 1:dim(cluster.table.live)[1]){
    first.pattern = as.character(cluster.table.live[b,3])
    if (first.pattern == cluster.table.live[b,4] && first.pattern == cluster.table.live[b,5] && first.pattern == cluster.table.live[b,6] && first.pattern == cluster.table.live[b,7] && first.pattern == cluster.table.live[b,8]){
        core.of.core = rbind(core.of.core, cluster.table.live[b,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
core.of.core.EEs.removed = c()
for (k in 1:dim(core.of.core)[1]){ #1, 2, 3, ... 247
    for (m in seq(3, 8, 1)){ #3, 4, ..., 8
        indicator <- isTRUE(core.of.core[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            core.of.core.EEs.removed = rbind(core.of.core.EEs.removed, core.of.core[k,])
            break
        }
    }
}
write.table(core.of.core.EEs.removed, file = "binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.the.same.across.infection.conditions.txt", quote=F, col.names=T, row.names = F) #22 genes


#4. What are the genes that behave the same way in the Gram-positive infection (M.luteus, E.fae.live)?
#Part 1. Pick out genes that have the same expression paths in both M.luteus and E.fae.live
cluster.table.positive = cluster.table[,c(1,2,3,6)]
gram.positive.only=c()
for (c in 1:dim(cluster.table.positive)[1]){
    first.pattern = as.character(cluster.table.positive[c,3]) #M.luteus pattern
    if (first.pattern == cluster.table.positive[c,4]){
        gram.positive.only = rbind(gram.positive.only, cluster.table.positive[c,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
gram.positive.only.EEs.removed = c()
for (k in 1:dim(gram.positive.only)[1]){ #1, 2, 3, ... 1295
    for (m in seq(3, 4, 1)){ #3, 4
        indicator <- isTRUE(gram.positive.only[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            gram.positive.only.EEs.removed = rbind(gram.positive.only.EEs.removed, gram.positive.only[k,])
            break
        }
    }
}
write.table(gram.positive.only.EEs.removed, file = "binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.the.same.in.gram.positive.infections.txt", quote=F, col.names=T, row.names = F) 
#326 genes



#5. What are the genes that behave the same way in the Gram-negative infection (E.coli, S.marcescens Type, P.rettgeri, Ecc15)?
#Part 1. Pick out genes that have the same expression paths in E.coli, S.marcescens Type, P.rettgeri, Ecc15
cluster.table.negative = cluster.table[,c(1,2,4:5,7:8)]; gram.negative.only=c()
for (d in 1:dim(cluster.table.negative)[1]){
    first.pattern = as.character(cluster.table.negative[d,3]) #E.coli pattern
    if (first.pattern == cluster.table.negative[d,4] && first.pattern == cluster.table.negative[d,5] && first.pattern == cluster.table.negative[d,6]){
        gram.negative.only = rbind(gram.negative.only, cluster.table.negative[d,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
gram.negative.only.EEs.removed = c()
for (k in 1:dim(gram.negative.only)[1]){ #1, 2, 3, ... 632
    for (m in seq(3, 6, 1)){ #3, 4, 5, 6
        indicator <- isTRUE(gram.negative.only[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            gram.negative.only.EEs.removed = rbind(gram.negative.only.EEs.removed, gram.negative.only[k,])
            break
        }
    }
}
write.table(gram.negative.only.EEs.removed, file = "binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.the.same.in.gram.negative.infections.txt", quote=F, col.names=T, row.names = F)


#6. What are the genes that behave the same way in the Gram-positive infection (M.luteus, E.fae.live) and in the Gram-negative infection (E.coli, S.marcescens Type, P.rettgeri, Ecc15), but behave differently between those Gram-type of infections? (ex. Up-EE-EE in Gram-positives and EE-EE-EE in Gram-negative)
#Part 1. Pick out genes that exist in both gram.positive.only and gram.negative.only
gram.positive.only$gene_id = factor(gram.positive.only$gene_id); gram.negative.only$gene_id = factor(gram.negative.only$gene_id)
gram.negative.only.ids = factor(gram.negative.only$gene_id); gram.positive.gram.negative.overlap = c()
for (c in 1:length(gram.negative.only.ids)){
    pattern = as.character(gram.negative.only.ids[c])
    if (pattern %in% gram.positive.only$gene_id == TRUE){
        positive = gram.positive.only[which(gram.positive.only$gene_id == pattern),]
        negative = gram.negative.only[which(gram.negative.only$gene_id == pattern),]
        together = cbind(positive, negative[,c(3:6)])
        gram.positive.gram.negative.overlap = rbind(gram.positive.gram.negative.overlap, together)
    }
}
#282 genes
#Part 2. Pick out genes that have different paths between Gram-types (This also removes genes that have EE-EE-EE in ALL conditions.)
gram.positive.gram.negative.overlap.true = c()
for (i in 1:dim(gram.positive.gram.negative.overlap)[1]){
    if (as.character(gram.positive.gram.negative.overlap[i,3]) != as.character(gram.positive.gram.negative.overlap[i,5])) { #If the gene patterns are not the same
        gram.positive.gram.negative.overlap.true = rbind(gram.positive.gram.negative.overlap.true, gram.positive.gram.negative.overlap[i,])
    }
}
write.table(gram.positive.gram.negative.overlap.true, file = "binning.the.clustering.results/edgeR_unchallenged_vs_infected_genes.that.behave.the.same.within.each.gram.infection.but.differently.between.gram.infections.txt", quote=F, col.names=T, row.names = F)



#7. Check the expression of individual genes
rnaseq.table = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T)
rnaseq.table[which(rnaseq.table$gene.id == "FBgn0002563"),][,3]

#Is there a gene that behaves differently depending on infection conditions?
cluster.table.live.na.omit #23 genes. No, there's no gene that has a different pattern in each infection condition.



########################################################################################################
#Part II. Previous time point infected  vs Present time point dataset
########################################################################################################

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")


#1. What are the genes that behave differently when infecting bacteria is live or heatkilled?
cluster.table = read.table("edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path_NS_removed.txt", header=T)
cluster.table.E.fae = cluster.table[,c(1,2,6,9)]; cluster.table.P.rett = cluster.table[,c(1,2,7,10)]
cluster.table.E.fae$E.fae.live.path = factor(cluster.table.E.fae$E.fae.live.path); cluster.table.E.fae$E.fae.heatkilled.path = factor(cluster.table.E.fae$E.fae.heatkilled.path)
cluster.table.P.rett$P.rett.live.path = factor(cluster.table.P.rett$P.rett.live.path); cluster.table.P.rett$P.rett.heatkilled.path = factor(cluster.table.P.rett$P.rett.heatkilled.path)

#E.faecalis live vs E.faecalis heatkilled. 744 genes have different patterns between conditions
E.faecalis.live.vs.heatkilled.infection.diff.genes = c()
for (k in 1:dim(cluster.table.E.fae)[1]){
    if (as.character(cluster.table.E.fae[k,3]) != as.character(cluster.table.E.fae[k,4])) { #If the gene patterns are not the same
        E.faecalis.live.vs.heatkilled.infection.diff.genes = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, cluster.table.E.fae[k,])
    }
}
write.table(E.faecalis.live.vs.heatkilled.infection.diff.genes, file="binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.differently.between.E.faecalis.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#P.rettgeri live vs P.rettgeri heatkilled. 1089 genes have different patterns between conditions
P.rettgeri.live.vs.heatkilled.infection.diff.genes = c()
for (n in 1:dim(cluster.table.P.rett)[1]){
    if (as.character(cluster.table.P.rett[n,3]) != as.character(cluster.table.P.rett[n,4])) { #If the gene patterns are not the same
        P.rettgeri.live.vs.heatkilled.infection.diff.genes = rbind(P.rettgeri.live.vs.heatkilled.infection.diff.genes, cluster.table.P.rett[n,])
    }
}
write.table(P.rettgeri.live.vs.heatkilled.infection.diff.genes, file="binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.differently.between.P.rettgeri.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)

#2. Are there common patterns of live vs heatkilled bacteria?
#Part 1. Pick genes that have changed between live and heatkilled in both infections.
colnames(E.faecalis.live.vs.heatkilled.infection.diff.genes)[3:4]= c("live","heatkilled")
colnames(P.rettgeri.live.vs.heatkilled.infection.diff.genes)[3:4]= c("live","heatkilled") #unify the column names
cluster.table.heat = rbind(E.faecalis.live.vs.heatkilled.infection.diff.genes, P.rettgeri.live.vs.heatkilled.infection.diff.genes)
dup_idx = duplicated(cluster.table.heat[,1]); cluster.table.heat.with.idx = cbind(cluster.table.heat, dup_idx) #pick genes that were different between live and heatkilled conditions in each infection and that had information in both infections
cluster.table.heat.common = cluster.table.heat.with.idx[which(cluster.table.heat.with.idx$dup_idx == "TRUE"),]; cluster.table.heat.common = c(as.character(cluster.table.heat.common[,1]))

cluster.table.heat.common.table = c()
for(a in 1:length(cluster.table.heat.common)){
    line = cluster.table[which(cluster.table[,1] == cluster.table.heat.common[a]),]
    cluster.table.heat.common.table = rbind(cluster.table.heat.common.table, line)
}
#Part 2. Find out which genes have the same pattern between live and heatkilled in both bacterial conditions.
cluster.table.heat.common.table = cluster.table.heat.common.table[,c(1,2,6,9,7,10)]
cluster.table.heat.same.table = c()
for (b in 1:dim(cluster.table.heat.common.table)[1]){
    if (as.character(cluster.table.heat.common.table[b,3]) == as.character(cluster.table.heat.common.table[b,5])
        && as.character(cluster.table.heat.common.table[b,4]) == as.character(cluster.table.heat.common.table[b,6])
    ) { #If the gene patterns are not the same
        cluster.table.heat.same.table = rbind(cluster.table.heat.same.table, cluster.table.heat.common.table[b,])
    }
}

write.table(cluster.table.heat.same.table, file = "binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_common_genes.that.behave.differently.between.live.and.heatkilled.txt", quote=F, col.names=T, row.names = F)
#Answer: 205 genes.


#3. What are the genes that behave the same way across infection conditions? ("core of core genes")
#Part 1. Pick out genes that have the same expression path across all conditions
cluster.table.live = cluster.table[,c(1:8)]; core.of.core=c()
for (b in 1:dim(cluster.table.live)[1]){
    first.pattern = as.character(cluster.table.live[b,3])
    if (first.pattern == cluster.table.live[b,4] && first.pattern == cluster.table.live[b,5] && first.pattern == cluster.table.live[b,6] && first.pattern == cluster.table.live[b,7] && first.pattern == cluster.table.live[b,8]){
        core.of.core = rbind(core.of.core, cluster.table.live[b,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
core.of.core.EEs.removed = c()
for (k in 1:dim(core.of.core)[1]){ #1, 2, 3, ... 123
    for (m in seq(3, 8, 1)){ #3, 4, ..., 8
        indicator <- isTRUE(core.of.core[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            core.of.core.EEs.removed = rbind(core.of.core.EEs.removed, core.of.core[k,])
            break
        }
    }
}
write.table(core.of.core.EEs.removed, file = "binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.across.infection.conditions.txt", quote=F, col.names=T, row.names = F) #20 genes


#4. What are the genes that behave the same way in the Gram-positive infection (M.luteus, E.fae.live)?
#Part 1. Pick out genes that have the same expression paths in both M.luteus and E.fae.live
cluster.table.positive = cluster.table[,c(1,2,3,6)]
gram.positive.only=c()
for (c in 1:dim(cluster.table.positive)[1]){
    first.pattern = as.character(cluster.table.positive[c,3]) #M.luteus pattern
    if (first.pattern == cluster.table.positive[c,4]){
        gram.positive.only = rbind(gram.positive.only, cluster.table.positive[c,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
gram.positive.only.EEs.removed = c()
for (k in 1:dim(gram.positive.only)[1]){ #1, 2, 3, ... 1051
    for (m in seq(3, 4, 1)){ #3, 4
        indicator <- isTRUE(gram.positive.only[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            gram.positive.only.EEs.removed = rbind(gram.positive.only.EEs.removed, gram.positive.only[k,])
            break
        }
    }
}
write.table(gram.positive.only.EEs.removed, file = "binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.in.gram.positive.infections.txt", quote=F, col.names=T, row.names = F) 
#267 genes



#5. What are the genes that behave the same way in the Gram-negative infection (E.coli, S.marcescens Type, P.rettgeri, Ecc15)?
#Part 1. Pick out genes that have the same expression paths in E.coli, S.marcescens Type, P.rettgeri, Ecc15
cluster.table.negative = cluster.table[,c(1,2,4:5,7:8)]; gram.negative.only=c()
for (d in 1:dim(cluster.table.negative)[1]){
    first.pattern = as.character(cluster.table.negative[d,3]) #E.coli pattern
    if (first.pattern == cluster.table.negative[d,4] && first.pattern == cluster.table.negative[d,5] && first.pattern == cluster.table.negative[d,6]){
        gram.negative.only = rbind(gram.negative.only, cluster.table.negative[d,])
    }
}
#Part 2. Remove those that are all EE-EE-EE in ALL conditions
gram.negative.only.EEs.removed = c()
for (k in 1:dim(gram.negative.only)[1]){ #1, 2, 3, ... 446
    for (m in seq(3, 6, 1)){ #3, 4, 5, 6
        indicator <- isTRUE(gram.negative.only[k,m] != "EE-EE-EE") #When the gene is up or down
        if (indicator == T){
            gram.negative.only.EEs.removed = rbind(gram.negative.only.EEs.removed, gram.negative.only[k,])
            break
        }
    }
}
write.table(gram.negative.only.EEs.removed, file = "binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.in.gram.negative.infections.txt", quote=F, col.names=T, row.names = F)
#44 genes


#6. What are the genes that behave the same way in the Gram-positive infection (M.luteus, E.fae.live) and in the Gram-negative infection (E.coli, S.marcescens Type, P.rettgeri, Ecc15), but behave differently between those Gram-type of infections? (ex. Up-EE-EE in Gram-positives and EE-EE-EE in Gram-negative)
#Part 1. Pick out genes that exist in both gram.positive.only and gram.negative.only
gram.positive.only$gene_id = factor(gram.positive.only$gene_id); gram.negative.only$gene_id = factor(gram.negative.only$gene_id)
gram.negative.only.ids = factor(gram.negative.only$gene_id); gram.positive.gram.negative.overlap = c()
for (c in 1:length(gram.negative.only.ids)){
    pattern = as.character(gram.negative.only.ids[c])
    if (pattern %in% gram.positive.only$gene_id == TRUE){
        positive = gram.positive.only[which(gram.positive.only$gene_id == pattern),]
        negative = gram.negative.only[which(gram.negative.only$gene_id == pattern),]
        together = cbind(positive, negative[,c(3:6)])
        gram.positive.gram.negative.overlap = rbind(gram.positive.gram.negative.overlap, together)
    }
}
#158 genes
#Part 2. Pick out genes that have different paths between Gram-types (This also removes genes that have EE-EE-EE in ALL conditions.)
gram.positive.gram.negative.overlap.true = c()
for (i in 1:dim(gram.positive.gram.negative.overlap)[1]){
    if (as.character(gram.positive.gram.negative.overlap[i,3]) != as.character(gram.positive.gram.negative.overlap[i,5])) { #If the gene patterns are not the same
        gram.positive.gram.negative.overlap.true = rbind(gram.positive.gram.negative.overlap.true, gram.positive.gram.negative.overlap[i,])
    }
}
write.table(gram.positive.gram.negative.overlap.true, file = "binning.the.clustering.results/edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.within.each.gram.infection.but.differently.between.gram.infections.txt", quote=F, col.names=T, row.names = F)



#7. Following the expression paths of core genes:
cluster.table.live = cluster.table[,c(1:8)]
core.gene.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
core.gene.list = core.gene.list[,1]

