#Overlap between infection, cleanprick, and heatkilled
#October 24-25, 2016 (Mon-Tue)
#Joo Hyun Im (ji72)

# Goal: Draw a venn diagram of genes that belong to:
# - infection: if significantly differentially expressed in any live infection in any time point compared to unchallenged,
# - sterile wound control: if significantly differentially expressed in any clean prick in any time point compared to unchallenged,
# - heatkilled: if significantly differentially expressed in any heatkilled infection in any time point compared to unchallenged

rm(list=ls(all=TRUE)) #delete any previous entry

#Input (significant in at least one condition -- whether it be live infection, clean prick, or heatkilled):
total_de_list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #2589 x 64 (using old data) **using the old data from 2015**


#1. Split the dataset based on conditions:
cleanprick = total_de_list[,c(3:8)]; M.luteus = total_de_list[,c(9:14)]; E.coli = total_de_list[,c(15:20)]; S.mar.type = total_de_list[,c(21:26)]; E.fae.live = total_de_list[,c(27:32)]; P.rett.live =total_de_list[,c(33:38)]; Ecc15 = total_de_list[,c(39:44)]; S.aureus = total_de_list[,c(45,46)]; P.sneebia = total_de_list[,c(47,48)]; S.mar.Db11 = total_de_list[,c(49,50)]; P.ento = total_de_list[,c(51,52)]; E.fae.heatkilled = total_de_list[,c(53:58)]; P.rett.heatkilled =total_de_list[,c(59:64)]

#list.of.conditions <- list(cleanprick, M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, S.aureus, P.sneebia, S.mar.Db11, P.ento, E.fae.heatkilled, P.rett.heatkilled)


#2. Draw the DEGs for infection
list.of.live.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live, Ecc15, S.aureus, P.sneebia, S.mar.Db11, P.ento)
deg.for.live.infection = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
for (i in 1:dim(total_de_list)[1]){ #Number of genes: 1,2,...2589
    for (j in 1:10){  #Number of conditions: 10 (M.luteus to P.entomophila)
        data.to.check = as.matrix(as.data.frame(list.of.live.conditions[j]))
        data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
        count = length(grep("Y", data.to.check.for.a.gene))
        count.accumulated = count.accumulated + count #the total number of Ys in a given gene (regardless of time points)
        if (count > 0 ) {
            specificity=specificity + 1 #the total number of conditions that has at least Ys
        }
    }
    specificity.total = rbind(specificity.total,specificity)
    if (j == 10 && count.accumulated >0){ #at least DEG in one live infection
        deg.for.live.infection = rbind(deg.for.live.infection, total_de_list[i,])
    }
    count.accumulated = 0; specificity = 0 
}
list.of.deg.for.live.infection = drop.terms(deg.for.live.infection[,1])
write.table(list.of.deg.for.live.infection, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.live.infection.txt", quote=F, row.names = F, col.names = F)


#3. Draw the genes for clean prick
deg.for.cleanprick = vector(); count=0
for (i in 1:dim(total_de_list)[1]){ #Number of genes: 1,2,...2589
    data.to.check = as.matrix(as.data.frame(cleanprick))
    data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
    count = length(grep("Y", data.to.check.for.a.gene)) #the total number of Ys in a given gene (regardless of time points)
    if (count > 0){ #at least DEG in one live infection
        deg.for.cleanprick = rbind(deg.for.cleanprick, total_de_list[i,])
    }
}
list.of.deg.for.cleanprick = drop.terms(deg.for.cleanprick[,1])
write.table(list.of.deg.for.cleanprick, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.cleanprick.txt", quote=F, row.names = F, col.names = F)


#4. Draw the genes for heatkilled
list.of.heatkilled.conditions <- list(E.fae.heatkilled, P.rett.heatkilled)
deg.for.heatkilled.infection = vector(); count=0; count.accumulated=0; specificity = 0; specificity.total = vector()
for (i in 1:dim(total_de_list)[1]){ #Number of genes: 1,2,...2589
    for (j in 1:2){  #Number of conditions: 2
        data.to.check = as.matrix(as.data.frame(list.of.heatkilled.conditions[j]))
        data.to.check.for.a.gene = t(as.matrix(data.to.check[i,]))
        count = length(grep("Y", data.to.check.for.a.gene))
        count.accumulated = count.accumulated + count #the total number of Ys in a given gene (regardless of time points)
        if (count > 0 ) {
            specificity=specificity + 1 #the total number of conditions that has at least Ys
        }
    }
    specificity.total = rbind(specificity.total,specificity)
    if (j == 2 && count.accumulated >0){ #at least DEG in one live infection
        deg.for.heatkilled.infection = rbind(deg.for.heatkilled.infection, total_de_list[i,])
    }
    count.accumulated = 0; specificity = 0 
}
list.of.deg.for.heatkilled.infection = drop.terms(deg.for.heatkilled.infection[,1])
write.table(list.of.deg.for.heatkilled.infection, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.all.degs.for.heatkilled.txt", quote=F, row.names = F, col.names = F)


#5. Find overlap
length(list.of.deg.for.live.infection)
length(list.of.deg.for.cleanprick)
length(list.of.deg.for.heatkilled.infection)

#i. overlap between all three conditions
overlap.between.live.and.cleanprick=c()
for (i in 1:length(list.of.deg.for.live.infection)){
    pattern = as.character(list.of.deg.for.live.infection[i])
    if (pattern %in% list.of.deg.for.cleanprick){
        overlap.between.live.and.cleanprick = c(overlap.between.live.and.cleanprick, pattern)
    }
}
length(overlap.between.live.and.cleanprick)
write.table(overlap.between.live.and.cleanprick, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/overlap.between.live.infections.and.cleanprick.txt", quote=F, row.names = F, col.names = F)


length(overlap.between.live.and.cleanprick) #112
overlap.between.live.and.cleanprick.and.heatkilled=c()
for (i in 1:length(overlap.between.live.and.cleanprick)){
    pattern = as.character(overlap.between.live.and.cleanprick[i])
    if (pattern %in% list.of.deg.for.heatkilled.infection){
        overlap.between.live.and.cleanprick.and.heatkilled = c(overlap.between.live.and.cleanprick.and.heatkilled, pattern)
    }
}
length(overlap.between.live.and.cleanprick.and.heatkilled) #90
write.table(overlap.between.live.and.cleanprick.and.heatkilled, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/overlap.between.live.and.cleanprick.and.heatkilled.txt", quote=F, row.names = F, col.names = F)


#ii. overlap between cleanprick and heatkilled
overlap.between.cleanprick.and.heatkilled=c()
for (i in 1:length(list.of.deg.for.cleanprick)){
    pattern = as.character(list.of.deg.for.cleanprick[i])
    if (pattern %in% list.of.deg.for.heatkilled.infection){
        overlap.between.cleanprick.and.heatkilled = c(overlap.between.cleanprick.and.heatkilled, pattern)
    }
}
length(overlap.between.cleanprick.and.heatkilled) #91
write.table(overlap.between.cleanprick.and.heatkilled, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/overlap.between.cleanprick.and.heatkilled.txt", quote=F, row.names = F, col.names = F)

#iii. overlap between live infection and heatkilled
overlap.between.live.and.heatkilled=c()
for (i in 1:length(list.of.deg.for.live.infection)){
    pattern = as.character(list.of.deg.for.live.infection[i])
    if (pattern %in% list.of.deg.for.heatkilled.infection){
        overlap.between.live.and.heatkilled = c(overlap.between.live.and.heatkilled, pattern)
    }
}
length(overlap.between.live.and.heatkilled) #483
write.table(overlap.between.live.and.heatkilled, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/overlap.between.live.and.heatkilled.txt", quote=F, row.names = F, col.names = F)


#6. What are the genes that are DE in live infection but are not in clean prick (sterile wound condition) nor heatkilled? Are the 7+ core genes a part of this live-infection-only genes or are they sporadic?

#essentially: list.of.deg.for.live.infection.only = list.of.deg.for.live.infection - overlap.between.live.and.cleanprick - overlap.between.live.and.heatkilled

# list.of.deg.for.live.infection = as.vector(list.of.deg.for.live.infection)
# #remove overlap.between.live.and.cleanprick from list.of.deg.for.live.infection
# list.of.deg.for.live.infection.only =c(); pattern = c(); table=c()
# for (i in 1:length(list.of.deg.for.live.infection)){
#     pattern = as.character(list.of.deg.for.live.infection[i])
#     if (pattern %in% overlap.between.live.and.cleanprick){
#     }
#     else if (pattern %in% overlap.between.live.and.heatkilled){
#     }
#     else{
#         table = total_de_list[which(total_de_list$gene_id == pattern),]
#         list.of.deg.for.live.infection.only = rbind(list.of.deg.for.live.infection.only, table)
#     }
# }
# length(list.of.deg.for.live.infection.only)
# write.table(list.of.deg.for.live.infection.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.live.infection.only.txt", quote=F, row.names = F, col.names = F)

seven.plus.core.genes.up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
seven.plus.core.genes.down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
seven.plus.core.genes = rbind(seven.plus.core.genes.up, seven.plus.core.genes.down) #252
seven.plus.core.genes = c(as.character(seven.plus.core.genes[,1]))

#Find overlap between core genes and genes that were DEGs only in live infection
list.of.deg.for.live.infection.only = list.of.deg.for.live.infection.only[,1]
overlap.between.live.only.and.core=c()
for (i in 1:length(list.of.deg.for.live.infection.only)){
    pattern = as.character(list.of.deg.for.live.infection.only[i])
    if (pattern %in% seven.plus.core.genes){
        overlap.between.live.only.and.core = c(overlap.between.live.only.and.core, pattern)
    }
}
length(overlap.between.live.only.and.core) #105 genes (updated on 1/17/17)


#7. What are the genes that are DE in heatkilled but are not in clean prick (sterile wound condition) nor live infection?

#essentially: list.of.deg.for.heatkilled.only = list.of.deg.for.heatkilled.infection - overlap.between.live.and.heatkilled - overlap.between.cleanprick.and.heatkilled

# list.of.deg.for.heatkilled.infection = as.vector(list.of.deg.for.heatkilled.infection)
# list.of.deg.for.heatkilled.infection.only =c(); pattern = c(); table=c()
# for (i in 1:length(list.of.deg.for.heatkilled.infection)){
#     pattern = as.character(list.of.deg.for.heatkilled.infection[i])
#     if (pattern %in% overlap.between.live.and.heatkilled){
#     }
#     else if (pattern %in% overlap.between.cleanprick.and.heatkilled){
#     }
#     else{
#         table = total_de_list[which(total_de_list$gene_id == pattern),]
#         list.of.deg.for.heatkilled.infection.only = rbind(list.of.deg.for.heatkilled.infection.only, table)
#     }
# }
# length(list.of.deg.for.heatkilled.infection.only)
# write.table(list.of.deg.for.heatkilled.infection.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.heatkilled.infection.only.txt", quote=F, row.names = F, col.names = F)
# 
# 
# #8. What are the genes that are DE in heatkilled but are not in clean prick (sterile wound condition) nor live infection?
# #essentially: list.of.deg.for.cleanprick.only = list.of.deg.for.cleanprick - overlap.between.live.and.cleanprick - overlap.between.cleanprick.and.heatkilled
# 
# list.of.deg.for.cleanprick = as.vector(list.of.deg.for.cleanprick)
# list.of.deg.for.cleanprick.only =c(); pattern = c(); table=c()
# for (i in 1:length(list.of.deg.for.cleanprick)){
#     pattern = as.character(list.of.deg.for.cleanprick[i])
#     if (pattern %in% overlap.between.live.and.cleanprick){
#     }
#     else if (pattern %in% overlap.between.cleanprick.and.heatkilled){
#     }
#     else{
#         table = total_de_list[which(total_de_list$gene_id == pattern),]
#         list.of.deg.for.cleanprick.only = rbind(list.of.deg.for.cleanprick.only, table)
#     }
# }
# write.table(list.of.deg.for.cleanprick.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.cleanprick.only.txt", quote=F, row.names = F, col.names = F)


#What if the overlap is between core live infection (not just live infection), sterile wound and heatkilled? (3/9/2017)
#Run everything as it is except: list.of.deg.for.live.infection

seven.plus.core.genes.up = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", header=T)
seven.plus.core.genes.down = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", header=T)
seven.plus.core.genes = rbind(seven.plus.core.genes.up, seven.plus.core.genes.down) #252
seven.plus.core.genes = c(as.character(seven.plus.core.genes[,1]))
list.of.deg.for.live.infection = seven.plus.core.genes

length(list.of.deg.for.live.infection) #252 core genes
length(list.of.deg.for.cleanprick) #114
length(list.of.deg.for.heatkilled.infection) #648

#Easier way to do it:
overlap.between.live.and.heatkilled = intersect(list.of.deg.for.live.infection,list.of.deg.for.heatkilled.infection)
length(overlap.between.live.and.heatkilled) #142
overlap.between.live.and.cleanprick= intersect(list.of.deg.for.live.infection,list.of.deg.for.cleanprick)
length(overlap.between.live.and.cleanprick) #83
overlap.between.cleanprick.and.heatkilled= intersect(list.of.deg.for.cleanprick,list.of.deg.for.heatkilled.infection)
length(overlap.between.cleanprick.and.heatkilled) #91

live.only = setdiff(list.of.deg.for.live.infection, overlap.between.live.and.heatkilled)
live.only = setdiff(live.only, overlap.between.live.and.cleanprick)
length(live.only)

cp.only = setdiff(list.of.deg.for.cleanprick, overlap.between.cleanprick.and.heatkilled)
cp.only = setdiff(cp.only, overlap.between.live.and.cleanprick)
length(cp.only)

hk.only = setdiff(list.of.deg.for.heatkilled.infection, overlap.between.live.and.heatkilled)
hk.only = setdiff(hk.only, overlap.between.cleanprick.and.heatkilled)
length(hk.only)

all.overlap = intersect(overlap.between.live.and.heatkilled,overlap.between.live.and.cleanprick)
core.overlapping.with.mamps.but.not.with.injury = setdiff(overlap.between.live.and.heatkilled, all.overlap)


write.table(live.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.live.only.when.core_genes_were_used.txt", quote=F, row.names = F, col.names = F)
write.table(cp.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.cleanprick.only.when.core_genes_were_used.txt", quote=F, row.names = F, col.names = F)
write.table(hk.only, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.for.heatkilled.only.when.core_genes_were_used.txt", quote=F, row.names = F, col.names = F)
write.table(all.overlap, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.overlapped.between.all.three.when.core_genes_were_used.txt", quote=F, row.names = F, col.names = F)
write.table(core.overlapping.with.mamps.but.not.with.injury, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/specific.comparisons/overlap_between_live_and_cleanprick_and_heatkilled/list.of.degs.overlapped.between.core.and.mamps.but.not.with.injury.when.core_genes_were_used.txt", quote=F, row.names = F, col.names = F)
