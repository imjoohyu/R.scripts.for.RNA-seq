#edgeR Commands for Mega RNA-seq for basic comparison
#Date: September 21, 2015 (updated on October 21, 2015, re-updated on November 8, 2015)
#ji72

rm(list=ls(all=TRUE))
library("edgeR")

#****Updated on November 8th, 2015****
#1. Get normalized counts from all conditions using all genes including the ones with super low counts.
#The purpose of this is to get a table of eligible genes with normalized counts so that we can feed it to the internal database.

#(i) Load the data
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples.txt")
dim(CountTable) #17558 x 105

#(ii) Get the pieces of normalization together
count.table <- cbind(CountTable[,c(2:35,37:70,72:105)]) #everything excluding 0hr unmolested
list.three = c(rep(c(  rep("UC",each=3), rep("Clean.prick",each=3), rep("M.luteus",each=3), rep("E.coli",each=3), 
                                            rep("S.mar.type",each=3), rep("E.faecalis",each=3), rep("P.rettgeri",each=3), rep("Ecc15",each=3), 
                                            "S.aureus", "P.sneebia", "S.mar.DB11","P.entomophila", rep("E.faecalis.heatkilled",each=3), 
                                            rep("P.rettgeri.heatkilled",each=3)),3))
time.three = c(rep( c(c("zero", "zero", "zero"), rep(c("twelve","thirty.six","five.half"), 7),rep("twelve", 4),rep(c("twelve","thirty.six","five.half"),2)),3) )

design <- data.frame(row.names = colnames(count.table), treatment = list.three, time = time.three, rep = c(rep(1:3, each = 34)) )
group <- factor(paste0(design$treatment, ".", design$time))

#(iii) Normalize the reads. Here, I decided to go for 4.4M *1.2= 5.28 count in at least expressed in 1/3 of the samples (34 samples). 4.4M is the average lib size
d = DGEList(counts = count.table, group=group); keep = rowSums(cpm(d) > 1.2) >=34; summary(keep);  d = d[keep, , keep.lib.sizes=FALSE]; 
d = calcNormFactors(d) #Estimate normalization factors
normalized.counts = cpm(d, normalized.lib.sizes=T) #Getting moderated log-counts-per-million from normalized counts from all conditions
normalized.counts.with.id = cbind(row.names(normalized.counts), normalized.counts)#[,c(1)]) #11911 genes after filtering ########
colnames(normalized.counts.with.id) = c("gene.id", "UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
                                        "M.luteus.12hr.Rep1","M.luteus.36hr.Rep1","M.luteus.5.5d.Rep1","E.coli.12hr.Rep1", "E.coli.36hr.Rep1","E.coli.5.5d.Rep1",
                                        "S.mar.type.12hr.Rep1","S.mar.type.36hr.Rep1","S.mar.type.5.5d.Rep1","E.fae.live.12hr.Rep1","E.fae.live.36hr.Rep1","E.fae.live.5.5d.Rep1",
                                        "P.rett.live.12hr.Rep1","P.rett.live.36hr.Rep1","P.rett.live.5.5d.Rep1","Ecc15.12hr.Rep1","Ecc15.36hr.Rep1","Ecc15.5.5d.Rep1",
                                        "S.aureus.12hr.Rep1","P.sneebia.12hr.Rep1","S.mar.Db11.12hr.Rep1","P.ento.12hr.Rep1","E.fae.heatkilled.12hr.Rep1","E.fae.heatkilled.36hr.Rep1","E.fae.heatkilled.5.5d.Rep1",
                                        "P.rett.heatkilled.12hr.Rep1","P.rett.heatkilled.36hr.Rep1","P.rett.heatkilled.5.5d.Rep1","UC.Rep4","UC.Rep5","UC.Rep6",
                                        "clean.prick.12hr.Rep2","clean.prick.36hr.Rep2","clean.prick.5.5d.Rep2",
                                        "M.luteus.12hr.Rep2","M.luteus.36hr.Rep2","M.luteus.5.5d.Rep2","E.coli.12hr.Rep2", "E.coli.36hr.Rep2","E.coli.5.5d.Rep2",
                                        "S.mar.type.12hr.Rep2","S.mar.type.36hr.Rep2","S.mar.type.5.5d.Rep2","E.fae.live.12hr.Rep2","E.fae.live.36hr.Rep2","E.fae.live.5.5d.Rep2",
                                        "P.rett.live.12hr.Rep2","P.rett.live.36hr.Rep2","P.rett.live.5.5d.Rep2","Ecc15.12hr.Rep2","Ecc15.36hr.Rep2","Ecc15.5.5d.Rep2",
                                        "S.aureus.12hr.Rep2","P.sneebia.12hr.Rep2","S.mar.Db11.12hr.Rep2","P.ento.12hr.Rep2","E.fae.heatkilled.12hr.Rep2","E.fae.heatkilled.36hr.Rep2","E.fae.heatkilled.5.5d.Rep2",
                                        "P.rett.heatkilled.12hr.Rep2","P.rett.heatkilled.36hr.Rep2","P.rett.heatkilled.5.5d.Rep2",
                                        "UC.Rep7","UC.Rep8","UC.Rep9","clean.prick.12hr.Rep3","clean.prick.36hr.Rep3","clean.prick.5.5d.Rep3",
                                        "M.luteus.12hr.Rep3","M.luteus.36hr.Rep3","M.luteus.5.5d.Rep3","E.coli.12hr.Rep3", "E.coli.36hr.Rep3","E.coli.5.5d.Rep3",
                                        "S.mar.type.12hr.Rep3","S.mar.type.36hr.Rep3","S.mar.type.5.5d.Rep3","E.fae.live.12hr.Rep3","E.fae.live.36hr.Rep3","E.fae.live.5.5d.Rep3",
                                        "P.rett.live.12hr.Rep3","P.rett.live.36hr.Rep3","P.rett.live.5.5d.Rep3","Ecc15.12hr.Rep3","Ecc15.36hr.Rep3","Ecc15.5.5d.Rep3",
                                        "S.aureus.12hr.Rep3","P.sneebia.12hr.Rep3","S.mar.Db11.12hr.Rep3","P.ento.12hr.Rep3","E.fae.heatkilled.12hr.Rep3","E.fae.heatkilled.36hr.Rep3","E.fae.heatkilled.5.5d.Rep3",
                                        "P.rett.heatkilled.12hr.Rep3","P.rett.heatkilled.36hr.Rep3","P.rett.heatkilled.5.5d.Rep3")

#(iv) Add gene names to the gene id #works
gene.id = normalized.counts.with.id[,1]; full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
output.mx <- matrix(NA, ncol=1, nrow=length(gene.id)[1]); colnames(output.mx) <-c("gene_name") #Create a matrix output that would hold gene names
for (i in 1:length(gene.id)[1]) {gene.to.compare = as.character(gene.id[i]); checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
output.mx[i,1] <- as.character(checking.against.full.list[1,10])    }
normalized.counts.with.name.and.id = cbind(normalized.counts.with.id[,1], output.mx, normalized.counts.with.id[,c(2:103)])
#write.table(normalized.counts.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_with_all_UCs.txt", quote=F, row.names=F)
write.table(normalized.counts.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_with_all_UCs_no_filter.txt", quote=F, row.names=F)


#(iv) Average the replicates
count = 1; normalized.counts.truncated = normalized.counts.with.name.and.id[,c(6:36,40:70,74:104)] #getting rid of gene ID and name and unchallenged samples
mean.table = matrix(NA, nrow=dim(normalized.counts.truncated)[1], ncol=31) #31 infection conditions
for (m in 1:dim(normalized.counts.truncated)[1]){
    for (n in 1:31){ #number of conditions is 31
        mean.value = mean(as.numeric(c(normalized.counts.truncated[m,count], normalized.counts.truncated[m,count+31], normalized.counts.truncated[m,count+62])))
        mean.table[m,n] = mean.value; count= count+1}    
        count=1}
normalized.counts.uc = normalized.counts.with.name.and.id[,c(3:5,37:39,71:73)];
mean.table2 = matrix(NA, nrow=dim(normalized.counts.uc)[1], ncol=1) #1 unchallenged condition
for (p in 1:dim(normalized.counts.uc)[1]){
    mean.value = mean(as.numeric(c(normalized.counts.uc[p,1], normalized.counts.uc[p,2], normalized.counts.uc[p,3],normalized.counts.uc[p,4], normalized.counts.uc[p,5], normalized.counts.uc[p,6],
                                   normalized.counts.uc[p,7], normalized.counts.uc[p,8], normalized.counts.uc[p,9])))
    mean.table2[p,1] = mean.value
}

mean.table.with.name.and.id = cbind(normalized.counts.with.name.and.id[,c(1:2)], mean.table2, mean.table)
colnames(mean.table.with.name.and.id) = c("gene.id","gene.name","UC","clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d","S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr","P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr","E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
#write.table(mean.table.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", quote=F, row.names=F, col.names=T)
write.table(mean.table.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs_no_filter.txt", quote=F, row.names=F, col.names=T)


#(v) Save the d$count as the raw data that has been filtered by cpm(d) > 1.2
raw.data.filtered = cbind(mean.table.with.name.and.id[,c(1:2)], d$counts) #11911 x 104 (everything except 0hr-unchallenged sample)
write.table(raw.data.filtered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", quote=F, row.names=F, col.names=T)

#******************
#2. Get DEGs from a basic comparison between an infection condition and the baseline (unchallenged) and filter out the genes accordingly
#Do the following job on CBSU
#(i) Load the data. Get the samples ready for DE analysis. I'm using 9 unchallenged samples (unchallenged 12h, 36hr, 5.5d x 3).
rm(list=ls(all=TRUE))
#setwd("/workdir/ji72")
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/")
library("edgeR")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T) #Use raw counts filtered by cpm earlier as an input
CountTable = CountTable[,c(3:104)] #Excluding the gene names and ids

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

#(ii) DE Analysis using edgeR: unchallenged vs infected
list.of.name.of.conditions <- c("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d",
                                "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
                                "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
list.of.conditions <- list(clean.prick.12hr, clean.prick.36hr, clean.prick.5.5d, M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d,
                           S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d, P.rett.live.12hr, P.rett.live.36hr,
                           P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr,
                           E.fae.heatkilled.12hr, E.fae.heatkilled.36hr, E.fae.heatkilled.5.5d, P.rett.heatkilled.12hr, P.rett.heatkilled.36hr, P.rett.heatkilled.5.5d)
de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1])

for (j in 1:length(list.of.conditions)){
    cat("edgeR - We are working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions[j])); count.table <- cbind(unchallenged.0hr, condition.count.table)
    if (j == 8){ #for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged",
                                                                                as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]))) } #removing 36hr rep2 sample
    else{
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged",
                                                                                as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]),as.character(list.of.name.of.conditions[j]))) }
    d = DGEList(counts = count.table, group=design.sp$condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est.commondisp = estimateCommonDisp(d); d.est.tagwisedisp = estimateTagwiseDisp(d.est.commondisp) #These two steps can be done in one step using estimateDisp() = equivalent
    
    #Plot a multidimensional scaling plot (MDS), a mean-variance relationship plot, BCV plot
    #png(file = paste("Mega_RNA-seq_edgeR_MDS_MV_BCV_plot_",list.of.name.of.conditions[j],".jpg", sep=""), width=1300, height=400)
    #par(mfrow = c(1,3)); plotMDS(d, labels=design.sp$row.names, col=c("darkgreen","blue")[factor(design.sp$condition)], cex=1)
    #plotMeanVar(d.est.tagwisedisp, show.tagwise.vars=TRUE, NBline=TRUE); plotBCV(d.est.tagwisedisp, cex=1); dev.off()
    
    #Test for differential expression ('classic' edgeR)
    de = exactTest(d.est.tagwisedisp, pair=c("Unchallenged", as.character(list.of.name.of.conditions[j])) ) #This is fine for 1:1 basic comparison
    tt = topTags(de, n=nrow(d.est.tagwisedisp)) #This command automatically sorts by the smallest FDR
    #tt = topTags(de, n=nrow(d.est.tagwisedisp), adjust.method="BH") #the same as above
    rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
    tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    #Plot a MA plot
    #png(file = paste("Mega_RNA-seq_edgeR_MA_plot_",list.of.name.of.conditions[j],".jpg", sep=""), width=500, height=500)
    #plotSmear(d.est.tagwisedisp, de.tags=deg, cex=1); dev.off()
    
    cat("edgeR - Finished! We are done working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
}


#(iii) Add gene name and id
de.of.all.conditions = de.of.all.conditions[,c(2:63)] #Get rid of the first column
print("Writing out fold change and FDR of mega RNA-seq edgeR results ...")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
total.DE.with.name.and.id = cbind(CountTable[,c(1:2)], de.of.all.conditions) 
colnames(total.DE.with.name.and.id) <- c("gene_id", "gene_name", "log2FC:clean.prick.12hr", "FDR:clean.prick.12hr","log2FC:clean.prick.36hr", "FDR:clean.prick.36hr",
                                         "log2FC:clean.prick.5.5d", "FDR:clean.prick.5.5d", "log2FC:M.luteus.12hr", "FDR:M.luteus.12hr", "log2FC:M.luteus.36hr", "FDR:M.luteus.36hr",
                                         "log2FC:M.luteus.5.5d", "FDR:M.luteus.5.5d", "log2FC:E.coli.12hr", "FDR:E.coli.12hr", "log2FC:E.coli.36hr", "FDR:E.coli.36hr",
                                         "log2FC:E.coli.5.5d", "FDR:E.coli.5.5d", "log2FC:S.mar.type.12hr", "FDR:S.mar.type.12hr", "log2FC:S.mar.type.36hr", "FDR:S.mar.type.36hr",
                                         "log2FC:S.mar.type.5.5d", "FDR:S.mar.type.5.5d", "log2FC:E.fae.live.12hr", "FDR:E.fae.live.12hr", "log2FC:E.fae.live.36hr", 
                                         "FDR:E.fae.live.36hr", "log2FC:E.fae.live.5.5d", "FDR:E.fae.live.5.5d", "log2FC:P.rett.live.12hr", "FDR:P.rett.live.12hr",
                                         "log2FC:P.rett.live.36hr", "FDR:P.rett.live.36hr", "log2FC:P.rett.live.5.5d", "FDR:P.rett.live.5.5d", "log2FC:Ecc15.12hr",
                                         "FDR:Ecc15.12hr","log2FC:Ecc15.36hr", "FDR:Ecc15.36hr", "log2FC:Ecc15.5.5d", "FDR:Ecc15.5.5d", "log2FC:S.aureus.12hr",
                                         "FDR:S.aureus.12hr", "log2FC:P.sneebia.12hr", "FDR:P.sneebia.12hr", "log2FC:S.mar.Db11.12hr", "FDR:S.mar.Db11.12hr",
                                         "log2FC:P.ento.12hr", "FDR:P.ento.12hr", "log2FC:E.fae.heatkilled.12hr", "FDR:E.fae.heatkilled.12hr", "log2FC:E.fae.heatkilled.36hr",
                                         "FDR:E.fae.heatkilled.36hr", "log2FC:E.fae.heatkilled.5.5d", "FDR:E.fae.heatkilled.5.5d", "log2FC:P.rett.heatkilled.12hr",
                                         "FDR:P.rett.heatkilled.12hr", "log2FC:P.rett.heatkilled.36hr", "FDR:P.rett.heatkilled.36hr", "log2FC:P.rett.heatkilled.5.5d",
                                         "FDR:P.rett.heatkilled.5.5d")
write.table(total.DE.with.name.and.id, file="edgeR_basic_comparison_all_genes_FC.txt", quote=F, row.names=F)  # -- for unchallenged


#(iv) Filtering the results based on some criteria
print("Filtering out the genes that have NA for p-val ...") #Filter out the genes that have NA for p-val
total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) #All of the cells have some non-NA values so there was no need for this command (unlike Deseq1/2) 

#(a) Only pick the genes that have significant p-vals in all of columns = genes that are changing expression significantly across ALL conditions.
print("Picking the genes that have sig p-vals in all of the columns ...")
sig.in.all.conditions <- vector()
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]) #row - 11911 genes
width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]) #column
indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11505
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.DE.with.name.and.id[k,m] > 0.05) #cutoff: FDR of 5%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break #non-significant p-values
    }
    if (m == width.of.table && indicator == F){
        #only saving the genes that have significant p-vals in ALL of the columns
        sig.in.all.conditions <- rbind(sig.in.all.conditions,total.DE.with.name.and.id[k,1:width.of.table])
    }
}
#Results: One gene (CG4757) is significantly DE in ALL of the conditions.

#(b) Change p-values to significant (Y) vs non-significant (N)
print("Converting the p-values to Y or N depending on their degree of significance ..."); total.with.names.sig.indicated <- total.DE.with.name.and.id
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]); width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.DE.with.name.and.id[k,m] > 0.05) #cutoff: FDR of 5%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) { #If the case is NOT significant,
            total.with.names.sig.indicated[k,m] = "N"
        }
        else {
            total.with.names.sig.indicated[k,m] = "Y"
        }
    }
}
write.table(total.with.names.sig.indicated, file="edgeR_basic_comparison_FDR_converted_to_Y-N.txt", quote=F, row.names=F)  # -- for unchallenged

#(c) Filter OUT the genes that have NON-significant p-val in all of columns
print("Filtering out the genes that have non-sig p-val in all of the columns ...")
filteredCountTable <- vector()
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1])
width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.with.names.sig.indicated[k,m] == "Y") #cutoff: FDR of 5%, Y = significant
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break #significant p-values
    }
    if (m == width.of.table && indicator == FALSE){
        filteredCountTable <- append(filteredCountTable, total.with.names.sig.indicated[k,1])   #only saving the gene ids that gave non-significant p-val in all of the columns
    }
}
cat(length(filteredCountTable), "is the number of genes that got filtered out")
#Results: filteredCountTable: there are 9322 genes whose FDR were non-significant in all conditions.

filtered.total.with.names <- total.with.names.sig.indicated[-c(filteredCountTable),]
dim(filtered.total.with.names) #2589 x 64
write.table(filtered.total.with.names, file="edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", quote=F, row.names=F)  # -- for unchallenged
filtered.total.with.names.with.pval <- total.DE.with.name.and.id[-c(filteredCountTable),]
dim(filtered.total.with.names.with.pval) #2589 x 64
write.table(filtered.total.with.names.with.pval, file="edgeR_basic_comparison_pval_at_least_one_sig.txt", quote=F, row.names=F)  # -- for unchallenged


#===========================
#October 2015

#1. Get normalized counts from all conditions using all genes including the ones with super low counts.
#The purpose of this is to get a table of eligible genes with normalized counts so that we can feed it to the internal database.

#(i) Load the data and average the unchallenged samples to 3 replicates
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples.txt")
Unchallenged.averaged = matrix(NA, nrow=dim(CountTable)[1], ncol=3) #Excluding unchallenged
for (i in 1:dim(CountTable)[1]){ #Averaging unchallenged samples
    Unchallenged.averaged[i,1] <- sum(CountTable[i,2],CountTable[i,3],CountTable[i,4])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,2] <- sum(CountTable[i,37],CountTable[i,38],CountTable[i,39])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
    Unchallenged.averaged[i,3] <- sum(CountTable[i,72],CountTable[i,73],CountTable[i,74])/3 #summing 12, 36, and 5.5 days unmolested readings and averaging by 3.
}
colnames(Unchallenged.averaged) <-c("UC.Rep1", "UC.Rep2", "UC.Rep3"); rownames(Unchallenged.averaged) <- row.names(CountTable)
Unchallenged.averaged <- round(Unchallenged.averaged, digits=0); count.table <- cbind(Unchallenged.averaged, CountTable[,c(5:35,40:70,75:105)]) #Adding averaged unchallenged

#(ii) Get the pieces of normalization together
list.three = c( c("UC", "UC", "UC"), rep(c( rep("Clean.prick",each=3), rep("M.luteus",each=3), rep("E.coli",each=3), 
                                            rep("S.mar.type",each=3), rep("E.faecalis",each=3), rep("P.rettgeri",each=3), rep("Ecc15",each=3), 
                                            "S.aureus", "P.sneebia", "S.mar.DB11","P.entomophila", rep("E.faecalis.heatkilled",each=3), 
                                            rep("P.rettgeri.heatkilled",each=3)),3))
time.three = c( c("zero", "zero", "zero"), rep( c(rep(c("twelve","thirty.six","five.half"), 7),rep("twelve", 4),rep(c("twelve","thirty.six","five.half"),2)),3) )
design <- data.frame(row.names = colnames(count.table), treatment = list.three, time = time.three, rep = c( c(1,2,3), rep(1:3, each = 31)) )
group <- factor(paste0(design$treatment, ".", design$time))

#(iii) Normalize the reads. Here, I decided to go for 4.4M *1.2= 5.28 count in at least expressed in 1/3 of the samples (35 samples). 4.4M is the average lib size
d = DGEList(counts = count.table, group=group); keep = rowSums(cpm(d) > 1.2) >=35; summary(keep);  d = d[keep, , keep.lib.sizes=FALSE]; d = calcNormFactors(d) #Estimate normalization factors
normalized.counts = cpm(d, normalized.lib.sizes=T) #Getting moderated log-counts-per-million from normalized counts from all conditions
normalized.counts.with.id = cbind(row.names(normalized.counts), normalized.counts)
colnames(normalized.counts.with.id) = c("gene.id", "UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1", "clean.prick.12hr.Rep2","clean.prick.12hr.Rep3","clean.prick.36hr.Rep1","clean.prick.36hr.Rep2","clean.prick.36hr.Rep3", "clean.prick.5.5d.Rep1", "clean.prick.5.5d.Rep2", "clean.prick.5.5d.Rep3","M.luteus.12hr.Rep1","M.luteus.12hr.Rep2","M.luteus.12hr.Rep3", "M.luteus.36hr.Rep1","M.luteus.36hr.Rep2","M.luteus.36hr.Rep3", "M.luteus.5.5d.Rep1", "M.luteus.5.5d.Rep2", "M.luteus.5.5d.Rep3","E.coli.12hr.Rep1","E.coli.12hr.Rep2","E.coli.12hr.Rep3", "E.coli.36hr.Rep1","E.coli.36hr.Rep2","E.coli.36hr.Rep3","E.coli.36hr.Rep1","E.coli.36hr.Rep2","E.coli.36hr.Rep3", "S.mar.type.12hr.Rep1","S.mar.type.12hr.Rep2","S.mar.type.12hr.Rep3","S.mar.type.36hr.Rep1","S.mar.type.36hr.Rep2","S.mar.type.36hr.Rep3","S.mar.type.5.5d.Rep1","S.mar.type.5.5d.Rep2","S.mar.type.5.5d.Rep3", "E.fae.live.12hr.Rep1","E.fae.live.12hr.Rep2","E.fae.live.12hr.Rep3", "E.fae.live.36hr.Rep1","E.fae.live.36hr.Rep2","E.fae.live.36hr.Rep3","E.fae.live.5.5d.Rep1","E.fae.live.5.5d.Rep2","E.fae.live.5.5d.Rep3","P.rett.live.12hr.Rep1","P.rett.live.12hr.Rep2","P.rett.live.12hr.Rep3","P.rett.live.36hr.Rep1","P.rett.live.36hr.Rep2","P.rett.live.36hr.Rep3","P.rett.live.5.5d.Rep1","P.rett.live.5.5d.Rep2","P.rett.live.5.5d.Rep3", "Ecc15.12hr.Rep1","Ecc15.12hr.Rep2","Ecc15.12hr.Rep3", "Ecc15.36hr.Rep1","Ecc15.36hr.Rep2","Ecc15.36hr.Rep3", "Ecc15.5.5d.Rep1","Ecc15.5.5d.Rep2","Ecc15.5.5d.Rep3", "S.aureus.12hr.Rep1","S.aureus.12hr.Rep2","S.aureus.12hr.Rep3","P.sneebia.12hr.Rep1","P.sneebia.12hr.Rep2","P.sneebia.12hr.Rep3", "S.mar.Db11.12hr.Rep1","S.mar.Db11.12hr.Rep2","S.mar.Db11.12hr.Rep3","P.ento.12hr.Rep1","P.ento.12hr.Rep2","P.ento.12hr.Rep3","E.fae.heatkilled.12hr.Rep1","E.fae.heatkilled.12hr.Rep2","E.fae.heatkilled.12hr.Rep3", "E.fae.heatkilled.36hr.Rep1","E.fae.heatkilled.36hr.Rep2","E.fae.heatkilled.36hr.Rep3", "E.fae.heatkilled.5.5d.Rep1","E.fae.heatkilled.5.5d.Rep2","E.fae.heatkilled.5.5d.Rep3", "P.rett.heatkilled.12hr.Rep1", "P.rett.heatkilled.12hr.Rep2","P.rett.heatkilled.12hr.Rep3", "P.rett.heatkilled.36hr.Rep1", "P.rett.heatkilled.36hr.Rep2","P.rett.heatkilled.36hr.Rep3","P.rett.heatkilled.5.5d.Rep1","P.rett.heatkilled.5.5d.Rep2","P.rett.heatkilled.5.5d.Rep3")

#(iv) Add gene names to the gene id #works
gene.id = normalized.counts.with.id[,1]; full.list <- read.table("/Users/JooHyun/Desktop/Drosophila_melanogaster.BDGP6.80.gene.list.txt",header=T)
output.mx <- matrix(NA, ncol=1, nrow=length(gene.id)[1]); colnames(output.mx) <-c("gene_name") #Create a matrix output that would hold gene names
for (i in 1:length(gene.id)[1]) {gene.to.compare = as.character(gene.id[i]); checking.against.full.list <- full.list[which(full.list$gene_id == gene.to.compare),]
output.mx[i,1] <- as.character(checking.against.full.list[1,10])    }
normalized.counts.with.name.and.id = cbind(normalized.counts.with.id[,1], output.mx, normalized.counts.with.id[,c(2:97)])
write.table(normalized.counts.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_Oct_2015/edgeR_normalized_counts_from_all_genes.txt", quote=F, row.names=F)

#(iv) Average the replicates
count = 1; normalized.counts.truncated = normalized.counts.with.name.and.id[,c(3:98)] #getting rid of gene ID and name
mean.table = matrix(NA, nrow=dim(normalized.counts.truncated)[1], ncol=32) #31 infection conditions + 1 averaged unchallenged = 32
for (m in 1:dim(normalized.counts.truncated)[1]){
    for (n in 1:32){ #number of conditions is 32
        mean.value = mean(as.numeric(c(normalized.counts.truncated[m,count], normalized.counts.truncated[m,count+1], normalized.counts.truncated[m,count+2])))
        mean.table[m,n] = mean.value; count= count+3}    
    count=1}
mean.table.with.name.and.id = cbind(normalized.counts.with.name.and.id[,c(1:2)], mean.table)
colnames(mean.table.with.name.and.id) = c("gene.id","gene.name","UC","clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d","S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr","P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr","E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
write.table(mean.table.with.name.and.id, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_Oct_2015/edgeR_normalized_counts_from_all_genes_averaged.txt", quote=F, row.names=F, col.names=T)

#(v) Save the d$count as the raw data that has been filtered by cpm(d) > 1.2
raw.data.filtered = cbind(mean.table.with.name.and.id[,c(1:2)], d$counts)
write.table(raw.data.filtered, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_Oct_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", quote=F, row.names=F, col.names=T)

#===========================

#2. Get DEGs from a basic comparison between an infection condition and the baseline (unchallenged) and filter out the genes accordingly
#Do the following job on CBSU
#(i) Load the data and average the unchallenged samples to 3 replicates. Get the samples ready for DE analysis
setwd("/workdir/ji72")
library("edgeR")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T) #Use raw counts filtered by cpm earlier as an input
CountTable = CountTable[,c(3:98)] #Excluding the gene names and ids

clean.prick.12hr = CountTable[,c(4,35,66)]; clean.prick.36hr = CountTable[,c(5,36,67)]; clean.prick.5.5d = CountTable[,c(6,37,68)]; M.luteus.12hr = CountTable[,c(7,38,69)]; M.luteus.36hr = CountTable[,c(8,39,70)]; M.luteus.5.5d = CountTable[,c(9,40,71)]; E.coli.12hr = CountTable[,c(10,41,72)]; E.coli.36hr = CountTable[,c(11,73)]
E.coli.5.5d = CountTable[,c(12,43,74)]; S.mar.type.12hr = CountTable[,c(13,44,75)]; S.mar.type.36hr = CountTable[,c(14,45,76)]; S.mar.type.5.5d = CountTable[,c(15,46,77)]; E.fae.live.12hr = CountTable[,c(16,47,78)]; E.fae.live.36hr = CountTable[,c(17,48,79)]; E.fae.live.5.5d = CountTable[,c(18,49,80)]; P.rett.live.12hr = CountTable[,c(19,50,81)]
P.rett.live.36hr = CountTable[,c(20,51,82)]; P.rett.live.5.5d = CountTable[,c(21,52,83)]; Ecc15.12hr = CountTable[,c(22,53,84)]; Ecc15.36hr = CountTable[,c(23,54,85)]; Ecc15.5.5d = CountTable[,c(24,55,86)]; S.aureus.12hr = CountTable[,c(25,56,87)]; P.sneebia.12hr = CountTable[,c(26,57,88)]; S.mar.Db11.12hr = CountTable[,c(27,58,89)]
P.ento.12hr = CountTable[,c(28,59,90)]; E.fae.heatkilled.12hr = CountTable[,c(29,60,91)]; E.fae.heatkilled.36hr = CountTable[,c(30,61,92)]; E.fae.heatkilled.5.5d = CountTable[,c(31,62,93)]; P.rett.heatkilled.12hr = CountTable[,c(32,63,94)]; P.rett.heatkilled.36hr = CountTable[,c(33,64,95)]; P.rett.heatkilled.5.5d = CountTable[,c(34,65,96)]

#list of conditions
list.of.name.of.conditions <- c("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d",
                                "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr",
                                "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr",
                                "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
list.of.conditions <- list(clean.prick.12hr, clean.prick.36hr, clean.prick.5.5d, M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d,
                           S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d, P.rett.live.12hr, P.rett.live.36hr,
                           P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr,
                           E.fae.heatkilled.12hr, E.fae.heatkilled.36hr, E.fae.heatkilled.5.5d, P.rett.heatkilled.12hr, P.rett.heatkilled.36hr, P.rett.heatkilled.5.5d)
de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1])

#(ii) DE Analysis using degeR
for (j in 1:length(list.of.conditions)){
    cat("edgeR - We are working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions[j])); count.table <- cbind(CountTable[,c(1:3)], condition.count.table)
    if (j == 8){ #for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Unchallenged","Unchallenged","Unchallenged",
                                                                                as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]))) } #removing 36hr rep2 sample
    else{
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Unchallenged","Unchallenged","Unchallenged",
                                                                            as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]),as.character(list.of.name.of.conditions[j]))) }
    d = DGEList(counts = count.table, group=design.sp$condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est.commondisp = estimateCommonDisp(d); d.est.tagwisedisp = estimateTagwiseDisp(d.est.commondisp) #These two steps can be done in one step using estimateDisp() = equivalent
    
    #Plot a multidimensional scaling plot (MDS), a mean-variance relationship plot, BCV plot
    png(file = paste("Mega_RNA-seq_edgeR_MDS_MV_BCV_plot_",list.of.name.of.conditions[j],".jpg", sep=""), width=1300, height=400)
    par(mfrow = c(1,3)); plotMDS(d, labels=design.sp$row.names, col=c("darkgreen","blue")[factor(design.sp$condition)], cex=1)
    plotMeanVar(d.est.tagwisedisp, show.tagwise.vars=TRUE, NBline=TRUE); plotBCV(d.est.tagwisedisp, cex=1); dev.off()
    
    #Test for differential expression ('classic' edgeR)
    de = exactTest(d.est.tagwisedisp, pair=c("Unchallenged", as.character(list.of.name.of.conditions[j])) ) #This is fine for 1:1 basic comparison
    tt = topTags(de, n=nrow(d.est.tagwisedisp)) #This command automatically sorts by the smallest FDR
    rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
    tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    #Plot a MA plot
    png(file = paste("Mega_RNA-seq_edgeR_MA_plot_",list.of.name.of.conditions[j],".jpg", sep=""), width=500, height=500)
    plotSmear(d.est.tagwisedisp, de.tags=deg, cex=1); dev.off()
    
    cat("edgeR - Finished! We are done working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
}

#(iii) Add gene name and id
de.of.all.conditions = de.of.all.conditions[,c(2:63)]
#Get rid of the first column
print("Writing out fold change and FDR of mega RNA-seq edgeR results ...")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
total.DE.with.name.and.id = cbind(CountTable[,c(1:2)], de.of.all.conditions) 
colnames(total.DE.with.name.and.id) <- c("gene_id", "gene_name", "log2FC:clean.prick.12hr", "FDR:clean.prick.12hr","log2FC:clean.prick.36hr", "FDR:clean.prick.36hr",
                                    "log2FC:clean.prick.5.5d", "FDR:clean.prick.5.5d", "log2FC:M.luteus.12hr", "FDR:M.luteus.12hr", "log2FC:M.luteus.36hr", "FDR:M.luteus.36hr",
                                    "log2FC:M.luteus.5.5d", "FDR:M.luteus.5.5d", "log2FC:E.coli.12hr", "FDR:E.coli.12hr", "log2FC:E.coli.36hr", "FDR:E.coli.36hr",
                                    "log2FC:E.coli.5.5d", "FDR:E.coli.5.5d", "log2FC:S.mar.type.12hr", "FDR:S.mar.type.12hr", "log2FC:S.mar.type.36hr", "FDR:S.mar.type.36hr",
                                    "log2FC:S.mar.type.5.5d", "FDR:S.mar.type.5.5d", "log2FC:E.fae.live.12hr", "FDR:E.fae.live.12hr", "log2FC:E.fae.live.36hr", 
                                    "FDR:E.fae.live.36hr", "log2FC:E.fae.live.5.5d", "FDR:E.fae.live.5.5d", "log2FC:P.rett.live.12hr", "FDR:P.rett.live.12hr",
                                    "log2FC:P.rett.live.36hr", "FDR:P.rett.live.36hr", "log2FC:P.rett.live.5.5d", "FDR:P.rett.live.5.5d", "log2FC:Ecc15.12hr",
                                    "FDR:Ecc15.12hr","log2FC:Ecc15.36hr", "FDR:Ecc15.36hr", "log2FC:Ecc15.5.5d", "FDR:Ecc15.5.5d", "log2FC:S.aureus.12hr",
                                    "FDR:S.aureus.12hr", "log2FC:P.sneebia.12hr", "FDR:P.sneebia.12hr", "log2FC:S.mar.Db11.12hr", "FDR:S.mar.Db11.12hr",
                                    "log2FC:P.ento.12hr", "FDR:P.ento.12hr", "log2FC:E.fae.heatkilled.12hr", "FDR:E.fae.heatkilled.12hr", "log2FC:E.fae.heatkilled.36hr",
                                    "FDR:E.fae.heatkilled.36hr", "log2FC:E.fae.heatkilled.5.5d", "FDR:E.fae.heatkilled.5.5d", "log2FC:P.rett.heatkilled.12hr",
                                    "FDR:P.rett.heatkilled.12hr", "log2FC:P.rett.heatkilled.36hr", "FDR:P.rett.heatkilled.36hr", "log2FC:P.rett.heatkilled.5.5d",
                                    "FDR:P.rett.heatkilled.5.5d")
write.table(total.DE.with.name.and.id, file="edgeR_basic_comparison_all_genes_FC.txt", quote=F, row.names=F)

#(iv) Filtering the results based on some criteria
print("Filtering out the genes that have NA for p-val ...") #Filter out the genes that have NA for p-val
total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) #All of the cells have some non-NA values so there was no need for this command (unlike Deseq1/2) 

#(a) Only pick the genes that have significant p-vals in all of columns = genes that are changing expression significantly across ALL conditions.
print("Picking the genes that have sig p-vals in all of the columns ...")
sig.in.all.conditions <- vector()
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]) #row - 11505 genes
width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]) #column
indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11505
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.DE.with.name.and.id[k,m] > 0.05) #cutoff: FDR of 5%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break #non-significant p-values
    }
    if (m == width.of.table && indicator == F){
        #only saving the genes that have significant p-vals in ALL of the columns
        sig.in.all.conditions <- rbind(sig.in.all.conditions,total.DE.with.name.and.id[k,1:width.of.table])
    }
}
#Results: No gene is significantly DE in ALL of the conditions.

#(b) Change p-values to significant (Y) vs non-significant (N)
print("Converting the p-values to Y or N depending on their degree of significance ..."); total.with.names.sig.indicated <- total.DE.with.name.and.id
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]); width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ... 11505
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.DE.with.name.and.id[k,m] > 0.05) #cutoff: FDR of 5%
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) { #If the case is NOT significant,
            total.with.names.sig.indicated[k,m] = "N"
        }
        else {
            total.with.names.sig.indicated[k,m] = "Y"
        }
    }
}
write.table(total.with.names.sig.indicated, file="edgeR_basic_comparison_FDR_converted_to_Y-N.txt", quote=F, row.names=F)

#(c) Filter OUT the genes that have NON-significant p-val in all of columns
print("Filtering out the genes that have non-sig p-val in all of the columns ...")
filteredCountTable <- vector()
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1])
width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11505
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 64
        #cat("The value for m is: ", m, "\n" )
        indicator <- isTRUE(total.with.names.sig.indicated[k,m] == "Y") #cutoff: FDR of 5%, Y = significant
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break #significant p-values
    }
    if (m == width.of.table && indicator == FALSE){
        filteredCountTable <- append(filteredCountTable, total.with.names.sig.indicated[k,1])   #only saving the gene ids that gave non-significant p-val in all of the columns
    }
}
cat(length(filteredCountTable), "is the number of genes that got filtered out")
#Results: filteredCountTable: there are 9677 genes whose FDR were non-significant in all conditions.

filtered.total.with.names <- total.with.names.sig.indicated[-c(filteredCountTable),]
dim(filtered.total.with.names) #2176 x 64
write.table(filtered.total.with.names, file="edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", quote=F, row.names=F)
filtered.total.with.names.with.pval <- total.DE.with.name.and.id[-c(filteredCountTable),]
dim(filtered.total.with.names.with.pval) #2176 x 64
write.table(filtered.total.with.names.with.pval, file="edgeR_basic_comparison_pval_at_least_one_sig.txt", quote=F, row.names=F)




#Leftover codes ...

# #Filter OUT the genes that have -1<log2FC<1 (2-fold difference) in all of columns
# print("Filtering out the genes that have -1<log2FC<1 in all of the columns ...")
# filteredCountTable.two <- vector(mode="integer", length=0)
# length.of.table <- as.numeric(dim(filtered.total.with.names)[1]) #row #from the command to get rid of non-significant genes 
# width.of.table <- as.numeric(dim(filtered.total.with.names)[2]) #column (64)
# indicator <- NULL
# for (k in 1:length.of.table){ #1, 2, 3, ...1887
#     #cat("The value for k is: ", k, "\n" )
#     for (m in seq(3, (width.of.table-1), 2)){ #3, 5, 7, ... 63
#         #cat("The value for m is: ", m, "\n" )
#         indicator <- isTRUE(
#             if (filtered.total.with.names[k,m] >0){
#             filtered.total.with.names[k,m] > 1 }
#             else{ 
#                 filtered.total.with.names[k,m] < -1}
#             ) #TRUE = what I want
#         #cat("indicator: ", indicator, "\n")
#         if (indicator == TRUE) break #FC is higher than 1 or lower than -1 
#     }
#     if (m == (width.of.table-1) && indicator == FALSE){
#         #print("I got to this point")
#         #print(filtered.total.with.names[k,1])
#         filteredCountTable.two <- append(filteredCountTable.two, filtered.total.with.names[k,1])   #only saving the gene ids that gave non-significant log2FC in all of the columns
#         #print(filteredCountTable.two)
#     }
# }

# #Results: filteredCountTable.two: there are 1311 genes whose log2FC were non-significant in all conditions.
# filtered.by.FDR.and.log2fc.with.names <- filtered.total.with.names[-filteredCountTable.two,]
# dim(filtered.by.FDR.and.log2fc.with.names) #576 x 64
# write.table(filtered.by.FDR.and.log2fc.with.names, file="Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_FCgene_ID-name_FDR_converted_to_Y-N_log2FC_filtered_at_least_one_sig.txt", quote=F, row.names=F)
# filtered.by.FDR.and.log2fc.with.names.with.pval <- filtered.total.with.names.with.pval[-c(filteredCountTable.two),]
# dim(filtered.by.FDR.and.log2fc.with.names.with.pval) #576 x 64
# write.table(filtered.by.FDR.and.log2fc.with.names.with.pval, file="Mega_RNA-seq_edgeR_basic_comparison_of_conditions_to_averaged_unchallenged_FCgene_ID-name_with_pval_log2FC_filtered_at_least_one_sig.txt", quote=F, row.names=F)
# 

#Find specific genes
# total.with.names.sig.indicated[total.with.names.sig.indicated[,2]=="Ald",] #Ald - example
# total.with.names.sig.indicated[total.with.names.sig.indicated[,2]=="CecB",] #Cecropin B
# total.with.names.sig.indicated[total.with.names.sig.indicated[,2]=="CecA",] #Cecropin A
# 
# #Look at which genes are significantly differentially expressed in ALL conditions
# sig.in.all.conditions

