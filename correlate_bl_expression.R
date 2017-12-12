#Correlate microbial load (BL) with expression level
#October 3rd, 2017
#Joo Hyun Im (ji72)

#David Schneider's (editor) comments:
#Is the fly responding to the number of microbes it faces or simply in a binary way in which it is infected or not? The microbe loads differ by 1000 fold; can this be correlated to the expression responses? Are there genes that only turn on only at high microbe loads? Are there genes that turn on at low microbe loads? What happens if you plot gene expression by microbe load for your favorite genes are there any interesting patterns that come out of this.

#Delete previous data
rm(list=ls(all=TRUE))


#Read in data
original_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
original_data_12h = original_data[c(1,2,9,15,21,27,33,39,45,47,49,51)] #Only choose 12h
colnames(original_data_12h) = c("gene_id", "gene_name", "Ml", "Ec", "SmType", "Ef", "Pr", "Ecc", "Sa", "Ps", "SmDb11", "PE")
#organize based on BL
original_data_12h = original_data_12h[,c(1:2,3,5,9,4,7,6,12,11,8,10)]


#Bacterial load at 12h (preliminary data from Katia)
bl_num = c(0, 3067.3, 8727.7, 11765, 39053.7, 66433.3, 117400, 715426.7, 909488.9, 1197333.3)
bl = data.frame(x=c(0,log10(bl_num[2:10])))
rownames(bl) = c("Ml", "SmType", "Sa", "Ec", "Pr", "Ef", "PE", "SmDb11", "Ecc",  "Ps")


#Calculate the correlation between microbial load and gene expression for each gene and sort the genes based on how highly correlated they are
cor_list = c()
for (i in 1:dim(original_data_12h)[1]){
    cor_list = rbind(cor_list, cor(t(original_data_12h[i,3:12]),bl$x))
}
original_data_12h_with_cor = cbind(original_data_12h, cor_list)
original_data_12h_with_cor_sorted = original_data_12h_with_cor[order(original_data_12h_with_cor$cor_list, decreasing = T),]

#When I looked at the first top 20, I didn’t see any immune genes.
head(original_data_12h_with_cor_sorted, 20)
xy = data.frame(cbind(t(original_data_12h[which(original_data_12h$gene_name == "Cyp6w1"),3:12]), bl)) #looking at a particular case
xy[order(xy[,1]),] #sorted by fold change

#When I specifically looked at immune genes, the correlation was not that great.
xy = data.frame(cbind(t(original_data_12h_with_cor_sorted[which(original_data_12h_with_cor_sorted$gene_name == "Dpt"),3:12]),bl)) #Dpt
colnames(xy) = c("fold_change", "bacterial_load"); xy[order(xy[,1]),] #sorted by fold change
xy = data.frame(cbind(t(original_data_12h_with_cor_sorted[which(original_data_12h_with_cor_sorted$gene_name == "Drs"),3:12]),bl)) #Drs
colnames(xy) = c("fold_change", "bacterial_load"); xy[order(xy[,1]),] #sorted by fold change

#Of all correlation values, highlight those with immune genes
plot(original_data_12h_with_cor_sorted$cor_list)
stripchart(original_data_12h_with_cor_sorted[which(original_data_12h_with_cor_sorted$gene_name %in% c("Dpt","Drs","PGRP-LC", "Rel", "PGRP-LB", "PGRP-SD", "AttA", "Mtk", "Spz", "Spirit", "PGRP-LA", "Spn88Eb", "TotA", "TotC", "TotM")),13], vertical =T, method="jitter",pch=10, col="maroon", add=TRUE)


#Are there genes that only turn on only at high microbe loads? Are there genes that turn on at low microbe loads?
#This is trickier to answer. Most of the top 20 genes that have a high correlation between bacterial load and gene expression patterns have less than 1 fold change upon infection. So I will exclude those genes that did not change expression very much first and reduce the number of genes to look at. The problem with this method is that it introduces biases.

#Keep if any of the conditions had > 0.5 or < -0.5 fold change
original_data_12h_with_cor_sorted_low_excluded = c()
for (m in 1:dim(original_data_12h_with_cor_sorted)[1]){
    for (n in 3:12){
        #I chose 0.5 as a cutoff because anecdotally 0.5 seems to be the lower bound of significance
        if (original_data_12h_with_cor_sorted[m,n] > 0.5 | original_data_12h_with_cor_sorted[m,n] < -0.5){
            original_data_12h_with_cor_sorted_low_excluded = rbind(original_data_12h_with_cor_sorted_low_excluded, original_data_12h_with_cor_sorted[m,])
            break
        }
    }
} #2174

#Another, more crude way
#original_data_12h_with_cor_sorted_low_excluded = original_data_12h_with_cor_sorted[which(original_data_12h_with_cor_sorted$PE > 0.8 | original_data_12h_with_cor_sorted$PE < -0.8),] #689

#If there are genes that turn on only at high microbe loads, the correlation between the gene and the microbe load should be highly positive. Likewise, if there are genes that turn on only at low microbe load, the correlation between the gene and the microbe load should be highly negative.

head(original_data_12h_with_cor_sorted_low_excluded,20)
tail(original_data_12h_with_cor_sorted_low_excluded,20)

plot(original_data_12h_with_cor_sorted_low_excluded$cor_list)
stripchart(original_data_12h_with_cor_sorted_low_excluded[which(original_data_12h_with_cor_sorted_low_excluded$gene_name %in% c("Dpt","Drs","PGRP-LC", "Rel", "PGRP-LB", "PGRP-SD", "AttA", "Mtk", "Spz", "Spirit", "PGRP-LA", "Spn88Eb", "TotA", "TotC", "TotM")),13], vertical =T, method="jitter",pch=10, col="maroon", add=TRUE)

#What are the genes whose expression levels are positively correlated with bacterial load? (cor>0.7) (GO terms)
original_data_12h_with_cor_sorted_low_excluded[which(original_data_12h_with_cor_sorted_low_excluded$cor_list >0.7),1] #24
write.table(original_data_12h_with_cor_sorted_low_excluded[which(original_data_12h_with_cor_sorted_low_excluded$cor_list >0.7),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Manuscript/RNAseq/Revision/correlate_bl_expression_highly_positive.txt", quote=F, row.names=F, col.names=T)

#What are the genes whose expression levels are negatively correlated with bacterial load? (cor>0.7) (GO terms)
original_data_12h_with_cor_sorted_low_excluded[which(original_data_12h_with_cor_sorted_low_excluded$cor_list < -0.7),1] #18
write.table(original_data_12h_with_cor_sorted_low_excluded[which(original_data_12h_with_cor_sorted_low_excluded$cor_list < -0.7),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Manuscript/RNAseq/Revision/correlate_bl_expression_highly_negative.txt", quote=F, row.names=F, col.names=T)



#Ml and Ecc stand out from the pack; Ml has no bacteria at 12h in all samples and Ecc has no bacteria at 12h in some samples.
#Previously I gave 0 for Ml's BL and 909488.9 for Ecc's BL. This might throw off some correlation.
#What if you remove Ml and Ecc from the list? What are the patterns then?

#Read in data
original_data_12h_subset = original_data_12h[c(1,2,4:10,12)] #Only choose 12h except Ml and Ecc15
bl_subset = data.frame(bl[c(2:8,10),]); rownames(bl_subset) = c("SmType", "Sa", "Ec", "Pr", "Ef", "PE", "SmDb11", "Ps")

#Calculate the correlation between microbial load and gene expression for each gene and sort the genes based on how highly correlated they are
cor_list_subset = c()
for (j in 1:dim(original_data_12h_subset)[1]){
    cor_list_subset = rbind(cor_list_subset, cor(t(original_data_12h_subset[j,3:10]), bl_subset))
}
original_data_12h_subset_with_cor = cbind(original_data_12h_subset, cor_list_subset)
colnames(original_data_12h_subset_with_cor)[11] = "cor_list_subset"
original_data_12h_subset_with_cor_sorted = original_data_12h_subset_with_cor[order(original_data_12h_subset_with_cor$cor_list_subset, decreasing = T),]

#The expression levels of immune genes do not seem to be highly correlated with bacterial load
head(original_data_12h_subset_with_cor_sorted, 20)
plot(original_data_12h_subset_with_cor_sorted$cor_list)
stripchart(original_data_12h_subset_with_cor_sorted[which(original_data_12h_subset_with_cor_sorted$gene_name %in% c("Dpt","Drs","PGRP-LC", "Rel", "PGRP-LB", "PGRP-SD", "AttA", "Mtk", "Spz", "Spirit", "PGRP-LA", "Spn88Eb", "TotA", "TotC", "TotM")),11], vertical =T, method="jitter",pch=10, col="maroon", add=TRUE)

#Remove low impact genes
#Keep if any of the conditions had > 0.5 or < -0.5 fold change
original_data_12h_subset_with_cor_sorted_low_excluded = c()
for (m in 1:dim(original_data_12h_subset_with_cor_sorted)[1]){
    for (n in 3:10){
        #I chose 0.5 as a cutoff because anecdotally 0.5 seems to be the lower bound of significance
        if (original_data_12h_subset_with_cor_sorted[m,n] > 0.5 | original_data_12h_subset_with_cor_sorted[m,n] < -0.5){
            original_data_12h_subset_with_cor_sorted_low_excluded = rbind(original_data_12h_subset_with_cor_sorted_low_excluded, original_data_12h_subset_with_cor_sorted[m,])
            break
        }
    }
} #2143

#more crude way
#original_data_12h_subset_with_cor_sorted_low_excluded = original_data_12h_subset_with_cor_sorted[which(original_data_12h_subset_with_cor_sorted$PE > 0.8 | original_data_12h_subset_with_cor_sorted$PE < -0.8),] #689

#If there are genes that turn on only at high microbe loads, the correlation between the gene and the microbe load should be highly positive. Likewise, if there are genes that turn on only at low microbe load, the correlation between the gene and the microbe load should be highly negative.

head(original_data_12h_subset_with_cor_sorted_low_excluded,20)
tail(original_data_12h_subset_with_cor_sorted_low_excluded,20)

plot(original_data_12h_subset_with_cor_sorted_low_excluded$cor_list_subset)
stripchart(original_data_12h_subset_with_cor_sorted_low_excluded[which(original_data_12h_subset_with_cor_sorted_low_excluded$gene_name %in% c("Dpt","Drs","PGRP-LC", "Rel", "PGRP-LB", "PGRP-SD", "AttA", "Mtk", "Spz", "Spirit", "PGRP-LA", "Spn88Eb", "TotA", "TotC", "TotM")),11], vertical =T, method="jitter",pch=10, col="maroon", add=TRUE)

#What are the genes whose expression levels are positively correlated with bacterial load? (cor>0.7) (GO terms)
original_data_12h_subset_with_cor_sorted_low_excluded[which(original_data_12h_subset_with_cor_sorted_low_excluded$cor_list_subset >0.7),1] #55
write.table(original_data_12h_subset_with_cor_sorted_low_excluded[which(original_data_12h_subset_with_cor_sorted_low_excluded$cor_list > 0.7),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Manuscript/RNAseq/Revision/correlate_bl_expression_excluding_ml_ecc_highly_positive.txt", quote=F, row.names=F, col.names=T)

#What are the genes whose expression levels are negatively correlated with bacterial load? (cor<-0.7) (GO terms)
original_data_12h_subset_with_cor_sorted_low_excluded[which(original_data_12h_subset_with_cor_sorted_low_excluded$cor_list_subset < -0.7),1] #25
write.table(original_data_12h_subset_with_cor_sorted_low_excluded[which(original_data_12h_subset_with_cor_sorted_low_excluded$cor_list < -0.7),], file="/Users/JooHyun/Dropbox/Cornell/Lab/Manuscript/RNAseq/Revision/correlate_bl_expression_excluding_ml_ecc_highly_negative.txt", quote=F, row.names=F, col.names=T)


#Draw a bacterial plot
bl_table = cbind(rownames(bl), bl)
bl_table$`rownames(bl)` = factor(bl_table$`rownames(bl)`, levels=c("Ml", "SmType", "Sa", "Ec", "Pr","Ef","PE","SmDb11","Ecc","Ps"))
plot(bl_table$x ~ bl_table$`rownames(bl)`, xlab="Bacteria", ylab="log2(Bacterial Load)")


#Added 12/5/2017
#My argument: 
#"However, a substantial number of genes continued to be differentially regulated at 36 h and 132 h post-inoculation (Fig 1C), presumably in part because the hosts continue to carry their bacterial infections at these later time points."

#Reviewer's comment:
#3. Line 189: I am not sure this was ever explicitly addressed – are differentially regulated genes at these later time points broadly correlated to whether or not individuals in the population, on average, suppressed their bacterial load or continued to have persistence?

#Previously, I've said:
#Some genes continue to be expressed in flies with persistent bacterial load while others go back to the basal expression level in flies that suppress bacterial load.
#However, we also know that some host responses continue to be induced regardless of whether bacteria persist or not into later time points. Ex. metabolism genes

#Problem with this answer:
#This is not precise enough.

#My interpretation of the question:
#a. Correlation between the number of DE genes at later time points and the conditions with varying degree of BL

rm(list=ls(all=TRUE))

#Number of genes
num_of_genes= read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/num.of.DEGs.final.txt", header=F, stringsAsFactors = F)

#Problem: We don't have BL at 36h and 132h from the actual experiments. We only have 12h, 24h, 48h, 6d from preliminary data.

#Time point 36h (with 24h BL)
#BL - 36h equivalent: Bacterial load at 24h (preliminary data from Katia), no 36h available
bl_num_24h = c(0, 1287.5, 4443.033, 18924.6, 1828820,  2144389)
bl_24h = data.frame(x=c(0,log10(bl_num_24h[2:6])))
rownames(bl_24h) = c("Ml",  "Ec", "SmType", "Ecc", "Ef", "Pr")

#Number of genes - 36h
num_of_genes_36h = num_of_genes[c(46:50,52,15:19,21),]
upregulated = c(161, 62, 139, 70, 137, 300)
cor(bl_24h, upregulated) #0.33
downregulated = c(66, 2, 41, 6, 357, 156)
cor(bl_24h, downregulated) #0.58
up_and_down_combined = upregulated + downregulated
cor(bl_24h, up_and_down_combined) #0.58144
bl_binary = c(0,1,1,1,1,1) #Categorical
cor(bl_binary, up_and_down_combined) #0.059

#Time point 36h (with 48h BL)
#BL - 36h equivalent: Bacterial load at 48h (preliminary data from Katia), no 36h available
bl_num_48h = c(0, 2834.2, 4578.367, 7081.4, 7301.667, 9911.366)
bl_48h = data.frame(x=c(0,log10(bl_num_48h[2:6])))
rownames(bl_48h) = c("Ml",  "Ecc", "Ef", "Ec", "Pr", "SmType")

#Number of genes - 48h
num_of_genes_48h = num_of_genes[c(46:50,52,15:19,21),]
upregulated = c(161, 62, 139, 70, 137, 300)
cor(bl_48h, upregulated) #-0.008
downregulated = c(66, 2, 41, 6, 357, 156)
cor(bl_48h, downregulated) #0.2022
up_and_down_combined = upregulated + downregulated
cor(bl_48h, up_and_down_combined) #0.144
bl_binary = c(0,1,1,1,1,1) #Categorical
cor(bl_binary, up_and_down_combined) #0.059

#Time point 132h:
#BL - 132h equivalent: Bacterial load at 6 day (preliminary data from Katia), no 132h available
bl_num_6d = c(0, 526.2333, 915.5, 3312.115, 3443.667, 9051)
bl_6d = data.frame(x=c(0,log10(bl_num_6d[2:6])))
rownames(bl_6d) = c("Ml", "Ecc", "Ec", "SmType", "Pr", "Ef")

#Number of genes - 132h
num_of_genes_132h = num_of_genes[c(55:59,61,24:28,30),]
upregulated = c(246, 111, 152, 214, 184, 168)
cor(bl_6d, upregulated) #-0.5022
downregulated = c(99, 20, 51, 59, 39, 33)
cor(bl_6d, downregulated) #-0.7901384
up_and_down_combined = upregulated + downregulated
cor(bl_6d, up_and_down_combined) #-0.6257
bl_binary = c(0,1,1,1,1,1) #Categorical
cor(bl_binary, up_and_down_combined) #-0.778



#Results:
#Binarily,
#At 36h, whether bacterial species sustain or not is not correlated with the number of DE genes (r^2=0.059).
#At 132h, whether bacterial species sustain or not is not strongly anti-correlated with the number of DE genes (r^2=-0.778)

#Quantitatively,
#At 36h, the higher bacterial load, the more DE genes there are, but the correlation is not that strong (r^2 = 0.58). 
#At 132h, the higher the bacterial load, the fewer DE genes there are, but the anti-correlation is not that strong (r^2 = -0.6257). 


#I will change the manuscript to:
#"However, a substantial number of genes continued to be differentially regulated at 36 h and 132 h post-inoculation (Fig 1C), presumably in part because the hosts continue to carry their bacterial infections at these later time points." "The number of differentially genes at these time points, though, was not strongly correlated to the amount of bacteria present in the populatin at that time (0.58 for 36h and -0.63 for 132h).

#In the response to reviewers, we could say:
#"M. luteus is the only bacterium that is not detected at 36h and 132h. Other bacteria are detected at those time points with varying degrees. The presence or the amount of bacteria present in flies is not strongly correlated with the number of DE genes for that bacterial condition at 36 h and is not strongly anti-correlated with the number of DE genes for that bacterial condition at 132h.


#3. Line 189: I am not sure this was ever explicitly addressed – are differentially regulated genes at these later time points broadly correlated to whether or not individuals in the population, on average, suppressed their bacterial load or continued to have persistence?



#b. Correlation between the degree of induction/repression at later time points and the conditions with varying degree of BL

#Delete previous data
rm(list=ls(all=TRUE))

#Bacterial load at 6 day (preliminary data from Katia), no 132h available
bl_num = c(0, 526.2333, 915.5, 3312.115, 3443.667, 9051)
bl = data.frame(x=c(0,log10(bl_num[2:6])))
rownames(bl) = c("Ml", "Ecc", "Ec","SmType", "Pr", "Ef")

#Read in data
original_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)
original_data_132h = original_data[c(1,2,13,19,25,31,37,43)] #Only choose 132h
colnames(original_data_132h) = c("gene_id", "gene_name", "Ml", "Ec", "SmType", "Ef", "Pr", "Ecc")
original_data_132h = original_data_132h[,c(1:2,3,8,4,5,7,6)] #organize based on BL

cor_list = c()
for (i in 1:dim(original_data_132h)[1]){
    cor_list = rbind(cor_list, cor(t(original_data_132h[i,3:8]),bl$x))
}
original_data_132h_with_cor = cbind(original_data_132h, cor_list)
original_data_132h_with_cor_sorted = original_data_132h_with_cor[order(original_data_132h_with_cor$cor_list, decreasing = T),] #2589

#Keep if any of the conditions had > 0.5 or < -0.5 fold change
original_data_132h_with_cor_sorted_low_excluded = c()
for (m in 1:dim(original_data_132h_with_cor_sorted)[1]){
    for (n in 3:8){
        #I chose 0.5 as a cutoff because anecdotally 0.5 seems to be the lower bound of significance
        if (original_data_132h_with_cor_sorted[m,n] > 0.5 | original_data_132h_with_cor_sorted[m,n] < -0.5){
            original_data_132h_with_cor_sorted_low_excluded = rbind(original_data_132h_with_cor_sorted_low_excluded, original_data_132h_with_cor_sorted[m,])
            break
        }
    }
} #1178

