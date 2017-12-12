#edgeR Commands for Mega RNA-seq for comparison to clean prick
#Date: December 2, 2015
#ji72

##MEASURING TRUE INFECTION RESPONSE##

#1. Get DEGs from a basic comparison between an infection condition and the clean prick condition and filter out the genes accordingly
#Do the following job on CBSU

#(i) Load the data and average the unchallenged samples to 3 replicates. Get the samples ready for DE analysis
rm(list=ls(all=TRUE))
setwd("/workdir/ji72")
library("edgeR")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T) #Use raw counts filtered by cpm earlier as an input
CountTable = CountTable[,c(3:104)] #Excluding the gene names and ids

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

#(ii) DE Analysis using edgeR: clean pricked vs infected
list.of.name.of.conditions.12hr <- c("M.luteus.12hr", "E.coli.12hr", "S.mar.type.12hr",  "E.fae.live.12hr", "P.rett.live.12hr", "Ecc15.12hr","S.aureus.12hr","P.sneebia.12hr","S.mar.Db11.12hr","P.ento.12hr","E.fae.heatkilled.12hr", "P.rett.heatkilled.12hr") #12
list.of.name.of.conditions.36hr <- c("M.luteus.36hr", "E.coli.36hr","S.mar.type.36hr","E.fae.live.36hr","P.rett.live.36hr", "Ecc15.36hr", "E.fae.heatkilled.36hr","P.rett.heatkilled.36hr") #8
list.of.name.of.conditions.5.5d <- c("M.luteus.5.5d", "E.coli.5.5d", "S.mar.type.5.5d", "E.fae.live.5.5d", "P.rett.live.5.5d","Ecc15.5.5d","E.fae.heatkilled.5.5d","P.rett.heatkilled.5.5d") #8
list.of.conditions.12hr <- list(M.luteus.12hr,E.coli.12hr, S.mar.type.12hr, E.fae.live.12hr,P.rett.live.12hr,Ecc15.12hr,S.aureus.12hr,P.sneebia.12hr,S.mar.Db11.12hr, P.ento.12hr,E.fae.heatkilled.12hr,P.rett.heatkilled.12hr)
list.of.conditions.36hr <- list(M.luteus.36hr,E.coli.36hr, S.mar.type.36hr,E.fae.live.36hr,P.rett.live.36hr,Ecc15.36hr,E.fae.heatkilled.36hr, P.rett.heatkilled.36hr)
list.of.conditions.5.5d <- list(M.luteus.5.5d,E.coli.5.5d,S.mar.type.5.5d,E.fae.live.5.5d, P.rett.live.5.5d,Ecc15.5.5d,E.fae.heatkilled.5.5d, P.rett.heatkilled.5.5d)
clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)]; clean.prick.5.5d = CountTable[,c(6,40,74)]

de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1]); for (j in 1:length(list.of.conditions.12hr)){ #change the name of 'list.of.conditions' to a specific time point
    cat("edgeR - We are working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.12hr[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions.12hr[j])); count.table <- cbind(clean.prick.12hr, condition.count.table)
    design.sp = data.frame(row.names = colnames(count.table), condition = c("Clean.pricked","Clean.pricked","Clean.pricked",as.character(list.of.name.of.conditions.12hr[j]), as.character(list.of.name.of.conditions.12hr[j]),as.character(list.of.name.of.conditions.12hr[j])))
    d = DGEList(counts = count.table, group=design.sp$condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est.commondisp = estimateCommonDisp(d); d.est.tagwisedisp = estimateTagwiseDisp(d.est.commondisp) #These two steps can be done in one step using estimateDisp() = equivalent
    
    #Plot a multidimensional scaling plot (MDS), a mean-variance relationship plot, BCV plot
    png(file = paste("Mega_RNA-seq_edgeR_MDS_MV_BCV_plot_",list.of.name.of.conditions.12hr[j],"_compared_to_cleanprick.jpg", sep=""), width=1300, height=400)
    par(mfrow = c(1,3)); plotMDS(d, labels=design.sp$row.names, col=c("darkgreen","blue")[factor(design.sp$condition)], cex=1)
    plotMeanVar(d.est.tagwisedisp, show.tagwise.vars=TRUE, NBline=TRUE); plotBCV(d.est.tagwisedisp, cex=1); dev.off()
    
    #Test for differential expression ('classic' edgeR)
    de = exactTest(d.est.tagwisedisp, pair=c("Clean.pricked", as.character(list.of.name.of.conditions.12hr[j])) ) #This is fine for 1:1 basic comparison
    tt = topTags(de, n=nrow(d.est.tagwisedisp)) #This command automatically sorts by the smallest FDR
    rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
    tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    #Plot a MA plot
    png(file = paste("Mega_RNA-seq_edgeR_MA_plot_",list.of.name.of.conditions.12hr[j],"_compared_to_cleanprick.jpg", sep=""), width=500, height=500)
    plotSmear(d.est.tagwisedisp, de.tags=deg, cex=1); dev.off()
    
    cat("edgeR - Finished! We are done working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.12hr[j]),"\n")
}; de.of.everything = de.of.all.conditions[,c(2:25)]
de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1]); for (j in 1:length(list.of.conditions.36hr)){ #change the name of 'list.of.conditions' to a specific time point
    cat("edgeR - We are working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.36hr[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions.36hr[j])); count.table <- cbind(clean.prick.36hr, condition.count.table)
    if (j == 2){ #for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Clean.pricked","Clean.pricked","Clean.pricked",as.character(list.of.name.of.conditions.36hr[j]), as.character(list.of.name.of.conditions.36hr[j]))) } #removing 36hr rep2 sample
    else{
        design.sp = data.frame(row.names = colnames(count.table), condition = c("Clean.pricked","Clean.pricked","Clean.pricked",as.character(list.of.name.of.conditions.36hr[j]), as.character(list.of.name.of.conditions.36hr[j]),as.character(list.of.name.of.conditions.36hr[j]))) }
    d = DGEList(counts = count.table, group=design.sp$condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est.commondisp = estimateCommonDisp(d); d.est.tagwisedisp = estimateTagwiseDisp(d.est.commondisp) #These two steps can be done in one step using estimateDisp() = equivalent
    
    #Plot a multidimensional scaling plot (MDS), a mean-variance relationship plot, BCV plot
    png(file = paste("Mega_RNA-seq_edgeR_MDS_MV_BCV_plot_",list.of.name.of.conditions.36hr[j],"_compared_to_cleanprick.jpg", sep=""), width=1300, height=400)
    par(mfrow = c(1,3)); plotMDS(d, labels=design.sp$row.names, col=c("darkgreen","blue")[factor(design.sp$condition)], cex=1)
    plotMeanVar(d.est.tagwisedisp, show.tagwise.vars=TRUE, NBline=TRUE); plotBCV(d.est.tagwisedisp, cex=1); dev.off()
    
    #Test for differential expression ('classic' edgeR)
    de = exactTest(d.est.tagwisedisp, pair=c("Clean.pricked", as.character(list.of.name.of.conditions.36hr[j])) ) #This is fine for 1:1 basic comparison
    tt = topTags(de, n=nrow(d.est.tagwisedisp)) #This command automatically sorts by the smallest FDR
    rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
    tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    #Plot a MA plot
    png(file = paste("Mega_RNA-seq_edgeR_MA_plot_",list.of.name.of.conditions.36hr[j],"_compared_to_cleanprick.jpg", sep=""), width=500, height=500)
    plotSmear(d.est.tagwisedisp, de.tags=deg, cex=1); dev.off()
    
    cat("edgeR - Finished! We are done working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.36hr[j]),"\n")
}; de.of.everything = cbind(de.of.everything, de.of.all.conditions[,c(2:17)])
de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1]); for (j in 1:length(list.of.conditions.5.5d)){ #change the name of 'list.of.conditions' to a specific time point
    cat("edgeR - We are working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.5.5d[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions.5.5d[j])); count.table <- cbind(clean.prick.5.5d, condition.count.table)
    design.sp = data.frame(row.names = colnames(count.table), condition = c("Clean.pricked","Clean.pricked","Clean.pricked",as.character(list.of.name.of.conditions.5.5d[j]), as.character(list.of.name.of.conditions.5.5d[j]),as.character(list.of.name.of.conditions.5.5d[j])))
    d = DGEList(counts = count.table, group=design.sp$condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est.commondisp = estimateCommonDisp(d); d.est.tagwisedisp = estimateTagwiseDisp(d.est.commondisp) #These two steps can be done in one step using estimateDisp() = equivalent
    
    #Plot a multidimensional scaling plot (MDS), a mean-variance relationship plot, BCV plot
    png(file = paste("Mega_RNA-seq_edgeR_MDS_MV_BCV_plot_",list.of.name.of.conditions.5.5d[j],"_compared_to_cleanprick.jpg", sep=""), width=1300, height=400)
    par(mfrow = c(1,3)); plotMDS(d, labels=design.sp$row.names, col=c("darkgreen","blue")[factor(design.sp$condition)], cex=1)
    plotMeanVar(d.est.tagwisedisp, show.tagwise.vars=TRUE, NBline=TRUE); plotBCV(d.est.tagwisedisp, cex=1); dev.off()
    
    #Test for differential expression ('classic' edgeR)
    de = exactTest(d.est.tagwisedisp, pair=c("Clean.pricked", as.character(list.of.name.of.conditions.5.5d[j])) ) #This is fine for 1:1 basic comparison
    tt = topTags(de, n=nrow(d.est.tagwisedisp)) #This command automatically sorts by the smallest FDR
    rn = rownames(tt$table); deg = rn[tt$table$FDR < 0.05] 
    tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    #Plot a MA plot
    png(file = paste("Mega_RNA-seq_edgeR_MA_plot_",list.of.name.of.conditions.5.5d[j],"_compared_to_cleanprick.jpg", sep=""), width=500, height=500)
    plotSmear(d.est.tagwisedisp, de.tags=deg, cex=1); dev.off()
    
    cat("edgeR - Finished! We are done working on the following comparison: Clean-pricked vs", as.character(list.of.name.of.conditions.5.5d[j]),"\n")
}; de.of.everything = cbind(de.of.everything, de.of.all.conditions[,c(2:17)])

de.of.all.conditions = de.of.everything[,c(1,2,25,26,41,42,3,4,27,28,43,44,5,6,29,30,45,46,7,8,31,32,47,48,9,10,33,34,49,50,11,12,35,36,51,52,13,14,15,16,17,18,19,20,21,22,37,38,53,54,23,24,39,40,55,56)] # in right order
print("Writing out fold change and FDR of mega RNA-seq edgeR results ...")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
total.DE.with.name.and.id = cbind(CountTable[,c(1:2)], de.of.all.conditions) 
colnames(total.DE.with.name.and.id) <- c("gene_id", "gene_name", "log2FC:M.luteus.12hr", "FDR:M.luteus.12hr", "log2FC:M.luteus.36hr", "FDR:M.luteus.36hr",
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
write.table(total.DE.with.name.and.id, file="edgeR_comparison_to_cleanprick_all_genes_FC.txt", quote=F, row.names=F)

#(iv) Filtering the results based on some criteria
print("Filtering out the genes that have NA for p-val ...") #Filter out the genes that have NA for p-val
total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) #All of the cells have some non-NA values so there was no need for this command (unlike Deseq1/2) 

#(a) Only pick the genes that have significant p-vals in all of columns = genes that are changing expression significantly across ALL conditions.
print("Picking the genes that have sig p-vals in all of the columns ...")
sig.in.all.conditions <- vector()
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]) #row - 11911 genes
width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]) #column
indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 58
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
for (k in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 58
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
write.table(total.with.names.sig.indicated, file="edgeR_comparison_to_cleanprick_FDR_converted_to_Y-N.txt", quote=F, row.names=F)

#(c) Filter OUT the genes that have NON-significant p-val in all of columns
print("Filtering out the genes that have non-sig p-val in all of the columns ...")
filteredCountTable <- vector()
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1])
width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ...11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 58
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
#Results: filteredCountTable: there are 10993 genes whose FDR were non-significant in all conditions.

filtered.total.with.names <- total.with.names.sig.indicated[-c(filteredCountTable),]
dim(filtered.total.with.names) # 918 x 58
write.table(filtered.total.with.names, file="edgeR_comparison_to_cleanprick_FDR_converted_to_Y-N_at_least_one_sig.txt", quote=F, row.names=F) 
filtered.total.with.names.with.pval <- total.DE.with.name.and.id[-c(filteredCountTable),]
dim(filtered.total.with.names.with.pval) #918 x 58
write.table(filtered.total.with.names.with.pval, file="edgeR_comparison_to_cleanprick_pval_at_least_one_sig.txt", quote=F, row.names=F)


#2. Pull out DE genes for a given condition
data.table = filtered.total.with.names
#data.table <- read.table("edgeR_comparison_to_cleanprick_FDR_converted_to_Y-N_at_least_one_sig.txt",header=T)
list.of.name.of.conditions <- list("M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d","S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", 
                                   "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr","P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", 
                                   "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr","E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", 
                                   "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
count=1
for (n in seq(4,58,2)){
    sig.genes.in.this.condition = subset(data.table, grepl("Y",data.table[,n]) ); print(dim(sig.genes.in.this.condition)[1])
    subset = sig.genes.in.this.condition[,c(1:2,n-1,n)]
    subset.sorted = subset[order(subset[,3], decreasing=T),]
    write.table(subset.sorted[which(subset.sorted[,3] > 0),], file=paste("list_of_upregulated_DE_genes_in_",list.of.name.of.conditions[count],".txt",sep=""), 
                quote=F, col.names= T, row.names=F)
    write.table(subset.sorted[which(subset.sorted[,3] < 0),], file=paste("list_of_downregulated_DE_genes_in_",list.of.name.of.conditions[count],".txt",sep=""), 
                quote=F, col.names= T, row.names=F)
    count = count +1
}

#3.Count the number of DEGs per infection-time condition using UNIX
# find /Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_comp_with_cleanprick_Dec_2015/list_of_DEGs -type f -exec wc -l {} + > num.of.DEGs.clean.prick.Dec.2015.txt
# Then edit the txt file with excel 
#   (ex. list_of_downregulated_DE_genes_in_E.fae.heatkilled.12hr.txt to down_in_E.fae.heatkilled.12hr)
#   Subtract 1 from the number of genes to account for the heading

#4. Make the graphs of the number of DEGs per infection-time condition
#read the file that has a number of DEGs per infection-time condition
num.degs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_comp_with_cleanprick_Dec_2015/list_of_DEGs/num.of.DEGs.clean.prick.Dec.2015.txt", header=F)
num.degs.down = num.degs[1:28,]; num.degs.up = num.degs[29:56,]

par(mfrow=c(2,3)) 
barplot(num.degs.up[1:12,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), ylim=c(0,700), xlab="12hr",cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.up[13:20,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.up[23:28,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot
barplot(num.degs.down[1:12,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), ylim=c(0,700), xlab="12hr", cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.down[13:20,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.down[23:28,2], col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot

#Just get the legend
barplot(c(1:12), col=c("pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), 
        legend = c("M.luteus","E.coli","S.mar.type","Ecc15","P.rett.live","P.rett.heatkilled",
                   "E.fae.live", "E.fae.heatkilled","S.aureus", "P.sneebia", "S.mar.Db11", "P.ento"))