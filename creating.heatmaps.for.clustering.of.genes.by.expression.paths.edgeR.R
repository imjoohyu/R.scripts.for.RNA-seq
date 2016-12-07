#Creating the heatmaps for clustering results based on edgeR expression path assignment
#December 7, 2016
#Joo Hyun Im (ji72)

#Remove previous entries
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/binning.the.clustering.results/")
library(gplots); library(RColorBrewer)

#A.Unchallenged vs Infected data
FC_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/edgeR_unchallenged_vs_infected_all_genes_FC.txt", header=T)

#A-1.Genes whose expression paths are conserved within Gram-positive live infections and within Gram-negative live infections but are different between Gram types
#M.luteus, E.faecalis, E.coli, S.marcescns Type, P.rettgeri, Ecc15
get_gene_list_and_data_A1 = function(gene_list, sample_list){
    list_of_genes = read.table(gene_list, header=T)
    list_of_genes_only = as.vector(list_of_genes$gene_id)
    FC_data_A1 = c()
    for (i in 1:length(list_of_genes_only)){
        FC_data_set = FC_data[which(FC_data$gene_id == list_of_genes_only[i]),sample_list]
        FC_data_A1 = rbind(FC_data_A1, FC_data_set)
    }
    row.names(FC_data_A1) = as.vector(list_of_genes$gene_name)
    return(as.matrix(FC_data_A1))
}
A1_data = get_gene_list_and_data_A1("edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.within.each.gram.infection.but.differently.between.gram.infections.txt", c(9,11,13,27,29,31,15,17,19,21,23,25,33,35,37,39,41,43))
colnames(A1_data) = c("M.luteus_12hr", "M.luteus_36hr", "M.luteus_5.5d","E.faecalis_12hr", "E.faecalis_36hr", "E.faecalis_5.5d","E.coli_12hr","E.coli_36hr","E.coli_5.5d","S.marcescns_12hr","S.marcescns_36hr","S.marcescns_5.5d","P.rettgeri_12hr","P.rettgeri_36hr","P.rettgeri_5.5d","Ecc15_12hr","Ecc15_36hr","Ecc15_5.5d")

#my_palette = colorRampPalette(brewer.pal(3,"Pastel1"))(n = 299)
my_palette = colorRampPalette(c("red","black","green"))(n = 299)
heatmap.2(A1_data, Colv=FALSE, density.info="none", dendrogram="row", trace="none", col=my_palette, margins = c(10, 5), scale=c("column"))

#A-2.Genes whose expression paths are conserved within live infections and within heatkilled infections but are different between the status of bacteria.
#E.fae.live, E.fae.heatkilled, P.rett.live, P.rett.heatkilled
get_gene_list_and_data_A2 = function(gene_list, sample_list){
    list_of_genes = read.table(gene_list, header=T)
    list_of_genes_only = as.vector(list_of_genes$gene_id)
    FC_data_A2 = c()
    for (i in 1:length(list_of_genes_only)){
        FC_data_set = FC_data[which(FC_data$gene_id == list_of_genes_only[i]),sample_list]
        FC_data_A2 = rbind(FC_data_A2, FC_data_set)
    }
    row.names(FC_data_A2) = as.vector(list_of_genes$gene_name)
    return(as.matrix(FC_data_A2))
}
A2_data = get_gene_list_and_data_A2("edgeR_prev_infected_vs_present_infected_common_genes.that.behave.differently.between.live.and.heatkilled.txt", c(27,29,31,53,55,57,33,35,37,59,61,63))
colnames(A2_data) = c("E.faecalis_12hr", "E.faecalis_36hr", "E.faecalis_5.5d","E.faecalis_hk__12hr", "E.faecalis_hk_36hr", "E.faecalis_hk_5.5d","P.rettgeri_12hr","P.rettgeri_36hr","P.rettgeri_5.5d","P.rettgeri_hk_12hr","P.rettgeri_hk_36hr","P.rettgeri_hk_5.5d")

#my_palette = colorRampPalette(brewer.pal(3,"Pastel1"))(n = 299)
my_palette = colorRampPalette(c("red","black","green"))(n = 299)
heatmap.2(A2_data, Colv=FALSE, density.info="none", dendrogram="row", trace="none", col=my_palette, margins = c(12, 12), srtRow = 45, scale=c("column"))

#A-3. Genes whose expression patterns are conserved across conditions (“core of core genes”).
#M.luteus, E.faecalis, E.coli, S.marcescns Type, P.rettgeri, Ecc15
get_gene_list_and_data_A3 = function(gene_list, sample_list){
    list_of_genes = read.table(gene_list, header=T)
    list_of_genes_only = as.vector(list_of_genes$gene_id)
    FC_data_A3 = c()
    for (i in 1:length(list_of_genes_only)){
        FC_data_set = FC_data[which(FC_data$gene_id == list_of_genes_only[i]),sample_list]
        FC_data_A3 = rbind(FC_data_A3, FC_data_set)
    }
    row.names(FC_data_A3) = as.vector(list_of_genes$gene_name)
    return(as.matrix(FC_data_A3))
}
A3_data = get_gene_list_and_data_A3("edgeR_prev_infected_vs_present_infected_genes.that.behave.the.same.across.infection.conditions.txt", c(9,11,13,27,29,31,15,17,19,21,23,25,33,35,37,39,41,43))
colnames(A3_data) = c("M.luteus_12hr", "M.luteus_36hr", "M.luteus_5.5d","E.faecalis_12hr", "E.faecalis_36hr", "E.faecalis_5.5d","E.coli_12hr","E.coli_36hr","E.coli_5.5d","S.marcescns_12hr","S.marcescns_36hr","S.marcescns_5.5d","P.rettgeri_12hr","P.rettgeri_36hr","P.rettgeri_5.5d","Ecc15_12hr","Ecc15_36hr","Ecc15_5.5d")

#my_palette = colorRampPalette(brewer.pal(3,"Pastel1"))(n = 299)
my_palette = colorRampPalette(c("red","black","green"))(n = 299)
heatmap.2(A3_data, Colv=FALSE, density.info="none", dendrogram="row", trace="none", col=my_palette, margins = c(12, 12), srtRow = 45, scale=c("column"))
