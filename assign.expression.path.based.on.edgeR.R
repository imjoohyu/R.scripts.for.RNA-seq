####################################################
#Assign expression paths to genes using the information from edgeR
#December 10, 2016
#Joo Hyun Im (ji72)
####################################################

#delete any previous input
rm(list=ls(all=TRUE))


#****************************************************************
#Part I: Make comparisons between 0hr-12hr, 0hr-36hr, 0hr-5.5d
#****************************************************************


#1. Get DEGs from a basic comparison between an infection condition and the baseline (unchallenged) and filter out the genes accordingly
####################################################
#Do the following job on CBSU --> On my local laptop (9/6/2016). CBSU gives me errors.

#Load the data. Get the samples ready for DE analysis. I'm using 9 unchallenged samples (unchallenged 12h, 36hr, 5.5d x 3).
#setwd("/workdir/ji72")
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

library("edgeR")
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T) #Use raw counts filtered by cpm earlier as an input
CountTable = CountTable[,c(3:104)] #Excluding the gene names and ids, so 11911 x 102

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)]; P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)]; P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

#2. unchallenged (0hr) vs infected (12hr, 36hr, 5.5d)
####################################################
#DE Analysis using edgeR: unchallenged (0hr) vs infected (12hr, 36hr, 5.5d)
list.of.name.of.conditions <- c("clean.prick.12hr", "clean.prick.36hr", "clean.prick.5.5d", "M.luteus.12hr", "M.luteus.36hr", "M.luteus.5.5d", "E.coli.12hr", "E.coli.36hr", "E.coli.5.5d", "S.mar.type.12hr", "S.mar.type.36hr", "S.mar.type.5.5d", "E.fae.live.12hr", "E.fae.live.36hr", "E.fae.live.5.5d", "P.rett.live.12hr", "P.rett.live.36hr", "P.rett.live.5.5d", "Ecc15.12hr", "Ecc15.36hr", "Ecc15.5.5d", "S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr", "E.fae.heatkilled.12hr", "E.fae.heatkilled.36hr", "E.fae.heatkilled.5.5d", "P.rett.heatkilled.12hr", "P.rett.heatkilled.36hr", "P.rett.heatkilled.5.5d")
list.of.conditions <- list(clean.prick.12hr, clean.prick.36hr, clean.prick.5.5d, M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d, E.coli.12hr, E.coli.36hr, E.coli.5.5d,S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d, E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d, P.rett.live.12hr, P.rett.live.36hr,P.rett.live.5.5d, Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d, S.aureus.12hr, P.sneebia.12hr, S.mar.Db11.12hr, P.ento.12hr, E.fae.heatkilled.12hr, E.fae.heatkilled.36hr, E.fae.heatkilled.5.5d, P.rett.heatkilled.12hr, P.rett.heatkilled.36hr, P.rett.heatkilled.5.5d)
de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1])

count.table=c()
for (j in 1:length(list.of.conditions)){
    cat("edgeR - We are working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
    condition.count.table <- as.matrix(as.data.frame(list.of.conditions[j])); count.table <- cbind(unchallenged.0hr, condition.count.table)
    cat(colnames(count.table),"\n")
    
    #for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr 
    if (j == 8){ 
        condition = c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged",as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j])) #removing 36hr rep2 sample
    } else { 
        condition = c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged",as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]), as.character(list.of.name.of.conditions[j]))
    }
    
    cat("condition: ", condition, "\n")
    design = model.matrix(~condition)
    d = DGEList(counts = count.table, group=condition) #Create a DGElist object
    d = calcNormFactors(d) #Estimate normalization factors
    d.est = estimateDisp(d, design)
    
    #Test for differential expression ('classic' edgeR)
    #de = exactTest(d.est, pair=c(as.character(condition[1]), as.character(condition[length(condition)])) ) #updated to this line (10/7/2016, 6:12pm)
    de = exactTest(d.est, pair=c("Unchallenged", as.character(list.of.name.of.conditions[j])) ) #This is fine for 1:1 basic comparison #old code as of 10/7/2016, 6:12pm #Back to this one as of 10/28/2016
    tt = topTags(de, n=nrow(d.est)) #This command automatically sorts by the smallest FDR
    rn = rownames(tt$table); tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
    
    #Save the output
    de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
    de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
    
    cat("edgeR - Finished! We are done working on the following comparison: Unchallenged vs", as.character(list.of.name.of.conditions[j]),"\n")
}

de.of.all.conditions = de.of.all.conditions[,c(2:63)] #Get rid of the first column
print("Writing out fold change and FDR of mega RNA-seq edgeR results ...")
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
total.DE.with.name.and.id = cbind(CountTable[,c(1:2)], de.of.all.conditions) 
colnames(total.DE.with.name.and.id) <- c("gene_id", "gene_name", "log2FC:clean.prick.12hr", "FDR:clean.prick.12hr","log2FC:clean.prick.36hr", "FDR:clean.prick.36hr","log2FC:clean.prick.5.5d", "FDR:clean.prick.5.5d", "log2FC:M.luteus.12hr", "FDR:M.luteus.12hr", "log2FC:M.luteus.36hr", "FDR:M.luteus.36hr","log2FC:M.luteus.5.5d", "FDR:M.luteus.5.5d", "log2FC:E.coli.12hr", "FDR:E.coli.12hr", "log2FC:E.coli.36hr", "FDR:E.coli.36hr","log2FC:E.coli.5.5d", "FDR:E.coli.5.5d", "log2FC:S.mar.type.12hr", "FDR:S.mar.type.12hr", "log2FC:S.mar.type.36hr", "FDR:S.mar.type.36hr","log2FC:S.mar.type.5.5d", "FDR:S.mar.type.5.5d", "log2FC:E.fae.live.12hr", "FDR:E.fae.live.12hr", "log2FC:E.fae.live.36hr", "FDR:E.fae.live.36hr", "log2FC:E.fae.live.5.5d", "FDR:E.fae.live.5.5d", "log2FC:P.rett.live.12hr", "FDR:P.rett.live.12hr","log2FC:P.rett.live.36hr", "FDR:P.rett.live.36hr", "log2FC:P.rett.live.5.5d", "FDR:P.rett.live.5.5d", "log2FC:Ecc15.12hr","FDR:Ecc15.12hr","log2FC:Ecc15.36hr", "FDR:Ecc15.36hr", "log2FC:Ecc15.5.5d", "FDR:Ecc15.5.5d", "log2FC:S.aureus.12hr","FDR:S.aureus.12hr", "log2FC:P.sneebia.12hr", "FDR:P.sneebia.12hr", "log2FC:S.mar.Db11.12hr", "FDR:S.mar.Db11.12hr","log2FC:P.ento.12hr", "FDR:P.ento.12hr", "log2FC:E.fae.heatkilled.12hr", "FDR:E.fae.heatkilled.12hr", "log2FC:E.fae.heatkilled.36hr", "FDR:E.fae.heatkilled.36hr", "log2FC:E.fae.heatkilled.5.5d", "FDR:E.fae.heatkilled.5.5d", "log2FC:P.rett.heatkilled.12hr","FDR:P.rett.heatkilled.12hr", "log2FC:P.rett.heatkilled.36hr", "FDR:P.rett.heatkilled.36hr", "log2FC:P.rett.heatkilled.5.5d","FDR:P.rett.heatkilled.5.5d")
write.table(total.DE.with.name.and.id, file="edgeR_unchallenged_vs_infected_all_genes_FC.txt", quote=F, row.names=F)  # -- for unchallenged

#Filter out the genes that have NA for p-val
print("Filtering out the genes that have NA for p-val ...") 
total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) #All of the cells have some non-NA values so there was no need for this command (unlike Deseq1/2) 

#Convert the p-values to Y or N depending on their significance
print("Converting the p-values to Y or N depending on their degree of significance ...")
total.with.names.sig.indicated <- total.DE.with.name.and.id
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
write.table(total.with.names.sig.indicated, file="edgeR_unchallenged_vs_infected_all_genes_FC_sig_changed_to_Y_or_N.txt", quote=F, row.names=F)


#Convert the FC information to Up, Down or Same depending on the direction of expression and the level of significance
print("Converting the FC information to Up, Down, or Same depending on their degree of significance ...")
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1]); width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2])
indicator <- NULL

for (i in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (s in seq(3, width.of.table, 2)){ #3, 5, 7, ... 63
        
        if (total.with.names.sig.indicated[i,s+1] == "Y"){ #if significant
            indicator <- isTRUE(total.DE.with.name.and.id[i,s] > 0) #indicator shows that the direction is positive
            #cat("indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the case is Up-DEG
                total.with.names.sig.indicated[i,s] = "Up"
            }
            else { #If the caseis Down-DEG
                total.with.names.sig.indicated[i,s] = "Down"
            }
        }
        else { #if not significant
            total.with.names.sig.indicated[i,s] = "EE"
        }
        
    }
}
write.table(total.with.names.sig.indicated, file="edgeR_unchallenged_vs_infected_all_genes_FC_converted_to_direction.txt", quote=F, row.names=F) 


#Put together the expression paths for the 6 live conditions and 2 heatkilled conditions
total.with.names.sig.indicated.path.only = total.with.names.sig.indicated[,c(9:44,53:64)] # 11911 x 48, excluding gene_id and gene_name
full.expression.path.table = matrix(NA, nrow = dim(total.with.names.sig.indicated.path.only)[1], ncol= 8)
#number.list = seq(3, dim(total.with.names.sig.indicated.path.only)[2], 2) #3, 5, 7, ..., 49

for (b in 1:dim(total.with.names.sig.indicated.path.only)[1]){ #for each gene
    count=1
    for (c in seq(1, (dim(total.with.names.sig.indicated.path.only)[2] -5),6)) { #for each condition of 6 live conditions and 2 heatkilled conditions 
        start = c; end = c+5    
        subset = total.with.names.sig.indicated.path.only[b,start:end]; subset
        path = paste(subset[1,1],"-",subset[1,3],"-",subset[1,5], sep="")
        full.expression.path.table[b, count] = path
        count=count+1
    }
}
full.expression.path.table = cbind(total.with.names.sig.indicated[,c(1:2)],full.expression.path.table)
colnames(full.expression.path.table) = c("gene_id","gene_name","M.luteus.path","E.coli.path","S.mar.type.path","E.fae.live.path","P.rett.live.path","Ecc15.path","E.fae.heatkilled.path","P.rett.heatkilled.path")
write.table(full.expression.path.table, file ="edgeR_unchallenged_vs_infected_all_genes_with_expression_path.txt", quote=F, row.names = F, col.names = T)

#Remove any genes that havent' changed significantly in any of the conditions.
full.expression.path.table.EEs.removed = c()
for (k in 1:dim(full.expression.path.table)[1]){ #1, 2, 3, ... 11911
    for (m in seq(3, 10, 1)){ #3, 4, ..., 10
        indicator <- isTRUE(full.expression.path.table[k,m] != "EE-EE-EE") #When the gene is up or down
        #cat("indicator: ", indicator, "\n")
        if (indicator == T){
            full.expression.path.table.EEs.removed = rbind(full.expression.path.table.EEs.removed, full.expression.path.table[k,])
            break
        }
    }
}
write.table(full.expression.path.table.EEs.removed, file ="edgeR_unchallenged_vs_infected_all_genes_with_expression_path_NS_removed.txt", quote=F, row.names = F, col.names = T)

#For a given condition, cluster genes based on patterns
M.luteus = full.expression.path.table[,c(1:3)]; E.coli = full.expression.path.table[,c(1:2,4)]
S.mar.type = full.expression.path.table[,c(1:2,5)]; E.fae.live = full.expression.path.table[,c(1:2,6)]
P.rett.live = full.expression.path.table[,c(1:2,7)]; Ecc15 = full.expression.path.table[,c(1:2,8)]
E.fae.heatkilled = full.expression.path.table[,c(1:2,9)]; P.rett.heatkilled = full.expression.path.table[,c(1:2,10)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

number.of.cluster=c()
for (g in 1:length(list.of.bacteria.name)){
    sample = data.frame(list.of.bacteria[g]) #Load in the expression path data
    number.of.cluster = c()
    sample[,3] = factor(sample[,3])
    pattern = levels(sample[,3])
    
    for (h in 1:length(pattern)){ #assign genes to clusters
    assign(paste0(list.of.bacteria.name[g],".cluster",h), sample[which(sample[,3] == pattern[g]),])
    }
    cat("Number of clusters for ", list.of.bacteria.name[g], " is ", length(pattern), "\n")
    number.of.cluster[g] = as.integer(length(pattern))
    
    sample.sorted = sample[order(sample[,3]),] #Sort the sample based on expression path and save them
    write.table(sample.sorted, file=paste("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_",list.of.bacteria.name[g],".txt",sep=""), quote=F, row.names = F, col.names = T )
}

#For a given condition, cluster genes based on patterns (Genes that were EE-EE-EE in ALL conditions were removed)
M.luteus = full.expression.path.table.EEs.removed[,c(1:3)]; E.coli = full.expression.path.table.EEs.removed[,c(1:2,4)]
S.mar.type = full.expression.path.table.EEs.removed[,c(1:2,5)]; E.fae.live = full.expression.path.table.EEs.removed[,c(1:2,6)]
P.rett.live = full.expression.path.table.EEs.removed[,c(1:2,7)]; Ecc15 = full.expression.path.table.EEs.removed[,c(1:2,8)]
E.fae.heatkilled = full.expression.path.table.EEs.removed[,c(1:2,9)]; P.rett.heatkilled = full.expression.path.table.EEs.removed[,c(1:2,10)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

number.of.cluster.EEs.removed=c()
for (g in 1:length(list.of.bacteria.name)){
    sample = data.frame(list.of.bacteria[g]) #Load in the expression path data
    number.of.cluster.EEs.removed = c()
    sample[,3] = factor(sample[,3])
    pattern = levels(sample[,3])
    
    for (h in 1:length(pattern)){ #assign genes to clusters
        assign(paste0(list.of.bacteria.name[g],".cluster",h), sample[which(sample[,3] == pattern[g]),])
    }
    cat("Number of clusters for ", list.of.bacteria.name[g], " is ", length(pattern), "\n")
    number.of.cluster.EEs.removed[g] = as.integer(length(pattern))
    
    sample.sorted = sample[order(sample[,3]),] #Sort the sample based on expression path and save them
    write.table(sample.sorted, file=paste("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_for_",list.of.bacteria.name[g],"_EEs_removed.txt",sep=""), quote=F, row.names = F, col.names = T )
}




#3. Plot the frequency of expression paths for each condition
####################################################
full.expression.path.table.EEs.removed = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path_NS_removed.txt", header=T) #2335 x 10
library(plyr)
M.luteus = full.expression.path.table.EEs.removed[,3]; M.luteus.counted = count(M.luteus); M.luteus.counted.sorted = M.luteus.counted[order(M.luteus.counted$x),]
E.coli = full.expression.path.table.EEs.removed[,4]; E.coli.counted = count(E.coli); E.coli.counted.sorted = E.coli.counted[order(E.coli.counted$x),]
S.mar.type = full.expression.path.table.EEs.removed[,5]; S.mar.type.counted = count(S.mar.type); S.mar.type.counted.sorted = S.mar.type.counted[order(S.mar.type.counted$x),]
Ecc15 = full.expression.path.table.EEs.removed[,8]; Ecc15.counted = count(Ecc15); Ecc15.counted.sorted = Ecc15.counted[order(Ecc15.counted$x),]
E.fae.live = full.expression.path.table.EEs.removed[,6]; E.fae.live.counted = count(E.fae.live); E.fae.live.counted.sorted = E.fae.live.counted[order(E.fae.live.counted$x),]
E.fae.heatkilled = full.expression.path.table.EEs.removed[,9]; E.fae.heatkilled.counted = count(E.fae.heatkilled); E.fae.heatkilled.counted.sorted = E.fae.heatkilled.counted[order(E.fae.heatkilled.counted$x),]
P.rett.live = full.expression.path.table.EEs.removed[,7]; P.rett.live.counted = count(P.rett.live); P.rett.live.counted.sorted = P.rett.live.counted[order(P.rett.live.counted$x),]
P.rett.heatkilled = full.expression.path.table.EEs.removed[,10]; P.rett.heatkilled.counted = count(P.rett.heatkilled); P.rett.heatkilled.counted.sorted = P.rett.heatkilled.counted[order(P.rett.heatkilled.counted$x),]

par(mfrow = c(2,4), mar=c(4,4,4,4), oma =c(1,5,0,0))
barplot(M.luteus.counted.sorted$freq, main = "M.luteus", horiz=T, names.arg=M.luteus.counted.sorted$x, las=1)
barplot(E.coli.counted.sorted$freq, main = "E.coli", horiz=T, names.arg=E.coli.counted.sorted$x, las=1)
barplot(S.mar.type.counted.sorted$freq, main = "S.marcescens Type", horiz=T, names.arg=S.mar.type.counted.sorted$x, las=1)
barplot(Ecc15.counted.sorted$freq, main = "Ecc15", horiz=T, names.arg=Ecc15.counted.sorted$x, las=1)
barplot(E.fae.live.counted.sorted$freq, main = "E.faecalis live", horiz=T, names.arg=E.fae.live.counted.sorted$x, las=1)
barplot(E.fae.heatkilled.counted.sorted$freq, main = "E.faecalis heatkilled", horiz=T, names.arg=E.fae.heatkilled.counted.sorted$x, las=1)
barplot(P.rett.live.counted.sorted$freq, main = "P.rettgeri live", horiz=T, names.arg=P.rett.live.counted.sorted$x, las=1)
barplot(P.rett.heatkilled.counted.sorted$freq, main = "P.rettgeri heatkilled", horiz=T, names.arg=P.rett.heatkilled.counted.sorted$x, las=1)

#Remove the genes that had EE-EE-EE in a given condition so that the barplot isn't skewed as much
M.luteus.counted.sorted = M.luteus.counted.sorted[which(M.luteus.counted.sorted$x != 'EE-EE-EE'),]
E.coli.counted.sorted = E.coli.counted.sorted[which(E.coli.counted.sorted$x != 'EE-EE-EE'),]
S.mar.type.counted.sorted = S.mar.type.counted.sorted[which(S.mar.type.counted.sorted$x != 'EE-EE-EE'),]
Ecc15.counted.sorted = Ecc15.counted.sorted[which(Ecc15.counted.sorted$x != 'EE-EE-EE'),]
E.fae.live.counted.sorted = E.fae.live.counted.sorted[which(E.fae.live.counted.sorted$x != 'EE-EE-EE'),]
E.fae.heatkilled.counted.sorted = E.fae.heatkilled.counted.sorted[which(E.fae.heatkilled.counted.sorted$x != 'EE-EE-EE'),]
P.rett.live.counted.sorted = P.rett.live.counted.sorted[which(P.rett.live.counted.sorted$x != 'EE-EE-EE'),]
P.rett.heatkilled.counted.sorted = P.rett.heatkilled.counted.sorted[which(P.rett.heatkilled.counted.sorted$x != 'EE-EE-EE'),]

par(mfrow = c(2,4), mar=c(4,4,4,4), oma =c(1,5,0,0))
barplot(M.luteus.counted.sorted$freq, main = "M.luteus", horiz=T, names.arg=M.luteus.counted.sorted$x, las=1, xlim=c(0,500))
barplot(E.coli.counted.sorted$freq, main = "E.coli", horiz=T, names.arg=E.coli.counted.sorted$x, las=1, xlim=c(0,500))
barplot(S.mar.type.counted.sorted$freq, main = "S.marcescens Type", horiz=T, names.arg=S.mar.type.counted.sorted$x, las=1, xlim=c(0,500))
barplot(Ecc15.counted.sorted$freq, main = "Ecc15", horiz=T, names.arg=Ecc15.counted.sorted$x, las=1, xlim=c(0,500))
barplot(E.fae.live.counted.sorted$freq, main = "E.faecalis live", horiz=T, names.arg=E.fae.live.counted.sorted$x, las=1, xlim=c(0,500))
barplot(E.fae.heatkilled.counted.sorted$freq, main = "E.faecalis heatkilled", horiz=T, names.arg=E.fae.heatkilled.counted.sorted$x, las=1, xlim=c(0,500))
barplot(P.rett.live.counted.sorted$freq, main = "P.rettgeri live", horiz=T, names.arg=P.rett.live.counted.sorted$x, las=1, xlim=c(0,500))
barplot(P.rett.heatkilled.counted.sorted$freq, main = "P.rettgeri heatkilled", horiz=T, names.arg=P.rett.heatkilled.counted.sorted$x, las=1, xlim=c(0,500))


#Create a list of expression paths
expression.list = read.table("edgeR_unchallenged_vs_infected_all_genes_with_expression_path.txt", header=T) #11911 x 10

#Pick the longest path
path.list = c()
for (i in c(3:10)){
    path.list = c(path.list, levels(expression.list[,i]))
}
path.list = unique(path.list)
path.list = data.frame(sort(path.list,decreasing = T)) #18 paths
write.table(path.list, file="list_of_expression_paths_from_unchallenged_vs_infected.txt", quote=F, row.names = T, col.names = F)







####################################################
#4. PCA of the expression path table, treating patterns as phenotypes: unchallenged (0hr) vs infected (12hr, 36hr, 5.5d)
full.expression.path.table.EEs.removed.path.only = full.expression.path.table.EEs.removed[,c(3:10)]
#Convert the paths into numbers (still has to be categorical variable)
path.list = read.table("list.of.expression.paths.txt", header=F) #this was compiled manually

for (i in 1:)


fit = factanal(full.expression.path.table.EEs.removed.path.only, 3, rotation="varimax")
print(fit, digits=2, cutoff=0.3, sort=T)
load = fit$loadings[,1:2]
plot(load, type="n"); text(load, labels=names(full.expression.path.table.EEs.removed.path.only), cex=1)






#****************************************************************
#Part II: Make comparisons between 0hr-12hr, 12hr-36hr, 36hr-5.5d
#****************************************************************
#1. Get DEGs from a basic comparison between the previous infection condition time point and the present infection condition time point and filter out the genes accordingly
#Do the following job on my local laptop (9/6/2016).


#delete any previous input
rm(list=ls(all=TRUE))

#Load the data. Get the samples ready for DE analysis.
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

library("edgeR")
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T) #Use raw counts filtered by cpm earlier as an input
CountTable = CountTable[,c(3:104)] #Excluding the gene names and ids

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)]; P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

####################################################
#2. previous time point infected vs next time point infected (0hr-12hr, 12hr-36hr, 36hr-5.5d)
#DE Analysis using edgeR

list.of.name.of.conditions <- c("Clean.prick", "M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")
Clean.prick.list = list(unchallenged.0hr,clean.prick.12hr,clean.prick.36hr,clean.prick.5.5d); M.luteus.list = list(unchallenged.0hr,M.luteus.12hr, M.luteus.36hr, M.luteus.5.5d); E.coli.list = list(unchallenged.0hr,E.coli.12hr, E.coli.36hr, E.coli.5.5d) 
S.mar.type.list = list(unchallenged.0hr,S.mar.type.12hr, S.mar.type.36hr, S.mar.type.5.5d); E.fae.live.list = list(unchallenged.0hr,E.fae.live.12hr, E.fae.live.36hr, E.fae.live.5.5d)
P.rett.live.list = list(unchallenged.0hr,P.rett.live.12hr, P.rett.live.36hr,P.rett.live.5.5d); Ecc15.list = list(unchallenged.0hr,Ecc15.12hr, Ecc15.36hr, Ecc15.5.5d)
E.fae.heatkilled.list = list(unchallenged.0hr,E.fae.heatkilled.12hr, E.fae.heatkilled.36hr, E.fae.heatkilled.5.5d); P.rett.heatkilled.list = list(unchallenged.0hr,P.rett.heatkilled.12hr, P.rett.heatkilled.36hr, P.rett.heatkilled.5.5d)

list.of.conditions <- list(Clean.prick.list, M.luteus.list, E.coli.list,S.mar.type.list,E.fae.live.list,P.rett.live.list,Ecc15.list,E.fae.heatkilled.list, P.rett.heatkilled.list)
time.points = c("0hr", "12hr", "36hr", "5.5d")

de.of.all.conditions <- matrix(NA, nrow=dim(CountTable)[1])
for (m in 1:length(list.of.name.of.conditions)){ #9 infection conditions
    for (n in 1:3){ #3 intervals between time points
        cat("edgeR - We are working on the following comparison: ", as.character(list.of.name.of.conditions[m])," ", time.points[n], " vs ", as.character(list.of.name.of.conditions[m]), " ", time.points[n+1], "\n")
        #cat("m: ", m, ", n: ",n)
        #Get the data
        count.table.1 = as.matrix(as.data.frame(list.of.conditions[[m]][[n]])) #time point 1
        count.table.2 = as.matrix(as.data.frame(list.of.conditions[[m]][[n+1]])) #time point 2
        count.table = cbind(count.table.1, count.table.2)
        cat(colnames(count.table))
        
        #Get the condition
        if (m == 3 && n == 2){ #12h3, 36h2, for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr when comparing 12hr and 36hr
            condition = c(paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."), paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."))
        } 
        else if (m == 3 && n == 3){ #36h2, 12h3, for E.coli 36hr -- getting rid of Rep2 of E.coli 36hr when comparing 36hr and 5.5d
            condition = c(paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."), paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."))
        }
        else if (n == 1){ #U9, 12h3
            condition = c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged","Unchallenged",paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."), paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."), paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."))
        }
        else if (n == 2 || n == 3){ #12h3, 36h3
            condition = c(paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."), paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."),paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep="."))
        }
        
        cat("condition: ",condition)
        design = model.matrix(~condition)
        d = DGEList(counts = count.table, group=condition) #Create a DGElist object
        d = calcNormFactors(d) #Estimate normalization factors
        d.est = estimateDisp(d, design)
        
        #Test for differential expression ('classic' edgeR)
        #The following code was adjusted on 10/28/2016 (Fri), 10:56am
        if (n == 1){ #when comparing UC(0hr) to 12hr
            de = exactTest(d.est, pair=c("Unchallenged",paste(as.character(list.of.name.of.conditions[m]),time.points[n+1],sep=".")))
        }
        
        else {
            de = exactTest(d.est, pair=c(as.character(condition[1]), as.character(condition[length(condition)])) ) #This is fine for 1:1 basic comparison
        }
        tt = topTags(de, n=nrow(d.est)) #This command automatically sorts by the smallest FDR
        rn = rownames(tt$table); tt.ordered.by.gene.id = tt$table[order(as.numeric(rownames(tt$table))),]
        
        #Save the output
        de.with.specific.info <- cbind(tt.ordered.by.gene.id[,1],tt.ordered.by.gene.id[,4]); colnames(de.with.specific.info) <- c("log2FC","FDR") #log2 Fold change and FDR
        de.of.all.conditions <- cbind(de.of.all.conditions, de.with.specific.info)
        
        cat("edgeR - Finished! We are done working on the following comparison: ", as.character(list.of.name.of.conditions[m])," ", time.points[n], " vs ", as.character(list.of.name.of.conditions[m])," ", time.points[n+1],"\n")
    }
}

de.of.all.conditions = de.of.all.conditions[,c(2:55)] #Get rid of the first column

print("Writing out fold change and FDR of mega RNA-seq edgeR results ...")
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/mega_RNA-seq_count_of_all_samples_filtered_cpm.txt", header=T)
total.DE.with.name.and.id = cbind(CountTable[,c(1:2)], de.of.all.conditions) 
colnames(total.DE.with.name.and.id) <- c("gene_id", "gene_name", "log2FC:clean.prick.12hr", "FDR:clean.prick.12hr","log2FC:clean.prick.36hr", "FDR:clean.prick.36hr",
                                         "log2FC:clean.prick.5.5d", "FDR:clean.prick.5.5d","log2FC:M.luteus.12hr", "FDR:M.luteus.12hr", "log2FC:M.luteus.36hr", "FDR:M.luteus.36hr",
                                         "log2FC:M.luteus.5.5d", "FDR:M.luteus.5.5d", "log2FC:E.coli.12hr", "FDR:E.coli.12hr", "log2FC:E.coli.36hr", "FDR:E.coli.36hr",
                                         "log2FC:E.coli.5.5d", "FDR:E.coli.5.5d", "log2FC:S.mar.type.12hr", "FDR:S.mar.type.12hr", "log2FC:S.mar.type.36hr", "FDR:S.mar.type.36hr",
                                         "log2FC:S.mar.type.5.5d", "FDR:S.mar.type.5.5d", "log2FC:E.fae.live.12hr", "FDR:E.fae.live.12hr", "log2FC:E.fae.live.36hr", 
                                         "FDR:E.fae.live.36hr", "log2FC:E.fae.live.5.5d", "FDR:E.fae.live.5.5d", "log2FC:P.rett.live.12hr", "FDR:P.rett.live.12hr",
                                         "log2FC:P.rett.live.36hr", "FDR:P.rett.live.36hr", "log2FC:P.rett.live.5.5d", "FDR:P.rett.live.5.5d", "log2FC:Ecc15.12hr",
                                         "FDR:Ecc15.12hr","log2FC:Ecc15.36hr", "FDR:Ecc15.36hr", "log2FC:Ecc15.5.5d", "FDR:Ecc15.5.5d","log2FC:E.fae.heatkilled.12hr", "FDR:E.fae.heatkilled.12hr", "log2FC:E.fae.heatkilled.36hr",
                                         "FDR:E.fae.heatkilled.36hr", "log2FC:E.fae.heatkilled.5.5d", "FDR:E.fae.heatkilled.5.5d", "log2FC:P.rett.heatkilled.12hr",
                                         "FDR:P.rett.heatkilled.12hr", "log2FC:P.rett.heatkilled.36hr", "FDR:P.rett.heatkilled.36hr", "log2FC:P.rett.heatkilled.5.5d",
                                         "FDR:P.rett.heatkilled.5.5d")
write.table(total.DE.with.name.and.id, file="edgeR_prev_infected_vs_present_infected_all_genes_FC.txt", quote=F, row.names=F)  # -- for unchallenged

#Filter out the genes that have NA for p-val
print("Filtering out the genes that have NA for p-val ...") 
total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) #All of the cells have some non-NA values so there was no need for this command (unlike Deseq1/2) 

#Convert the p-values to Y or N depending on their significance
print("Converting the p-values to Y or N depending on their degree of significance ...")
total.with.names.sig.indicated <- total.DE.with.name.and.id
length.of.table <- as.numeric(dim(total.DE.with.name.and.id)[1]); width.of.table <- as.numeric(dim(total.DE.with.name.and.id)[2]); indicator <- NULL
for (k in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (m in seq(4, width.of.table, 2)){ #4, 6, 8, ... 56
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
write.table(total.with.names.sig.indicated, file="edgeR_prev_infected_vs_present_infected_all_genes_FC_sig_changed_to_Y_or_N.txt", quote=F, row.names=F)


#Convert the FC information to Up, Down or Same depending on the direction of expression and the level of significance
print("Converting the FC information to Up, Down, or Same depending on their degree of significance ...")
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1]); width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2])
indicator <- NULL

for (i in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (s in seq(3, width.of.table, 2)){ #3, 5, 7, ... 56
        
        if (total.with.names.sig.indicated[i,s+1] == "Y"){ #if significant
            indicator <- isTRUE(total.DE.with.name.and.id[i,s] > 0) #indicator shows that the direction is positive
            #cat("indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the case is Up-DEG
                total.with.names.sig.indicated[i,s] = "Up"
            }
            else { #If the caseis Down-DEG
                total.with.names.sig.indicated[i,s] = "Down"
            }
        }
        else { #if not significant
            total.with.names.sig.indicated[i,s] = "EE"
        }
        
    }
}
write.table(total.with.names.sig.indicated, file="edgeR_prev_infected_vs_present_infected_all_genes_FC_converted_to_direction.txt", quote=F, row.names=F) 


#Put together the expression paths for the 6 live conditions and 2 heatkilled conditions
total.with.names.sig.indicated.path.only = total.with.names.sig.indicated[,c(9:56)] # 11911 x 56, excluding gene_id and gene_name and clean prick samples
full.expression.path.table = matrix(NA, nrow = dim(total.with.names.sig.indicated.path.only)[1], ncol= 8)

for (b in 1:dim(total.with.names.sig.indicated.path.only)[1]){ #for each gene
    count=1
    for (c in seq(1, (dim(total.with.names.sig.indicated.path.only)[2] -5),6)) { #for each condition of 6 live conditions and 2 heatkilled conditions 
        start = c; end = c+5    
        subset = total.with.names.sig.indicated.path.only[b,start:end]; subset
        path = paste(subset[1,1],"-",subset[1,3],"-",subset[1,5], sep="")
        full.expression.path.table[b, count] = path
        count=count+1
    }
}
full.expression.path.table = cbind(total.with.names.sig.indicated[,c(1:2)],full.expression.path.table)
colnames(full.expression.path.table) = c("gene_id","gene_name","M.luteus.path","E.coli.path","S.mar.type.path","E.fae.live.path","P.rett.live.path","Ecc15.path","E.fae.heatkilled.path","P.rett.heatkilled.path")
write.table(full.expression.path.table, file ="edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path.txt", quote=F, row.names = F, col.names = T)

#Remove any genes that havent' changed significantly. 1800 genes were left after filtering
full.expression.path.table.EEs.removed = c()
for (k in 1:dim(full.expression.path.table)[1]){ #1, 2, 3, ... 11911
    for (m in seq(3, 10, 1)){ #3, 4, ..., 10
        indicator <- isTRUE(full.expression.path.table[k,m] != "EE-EE-EE") #When the gene is up or down
        #cat("indicator: ", indicator, "\n")
        if (indicator == T){
            full.expression.path.table.EEs.removed = rbind(full.expression.path.table.EEs.removed, full.expression.path.table[k,])
            break
        }
    }
}
write.table(full.expression.path.table.EEs.removed, file ="edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path_NS_removed.txt", quote=F, row.names = F, col.names = T)

#For a given condition, cluster genes based on patterns
M.luteus = full.expression.path.table[,c(1:3)]; E.coli = full.expression.path.table[,c(1:2,4)]
S.mar.type = full.expression.path.table[,c(1:2,5)]; E.fae.live = full.expression.path.table[,c(1:2,6)]
P.rett.live = full.expression.path.table[,c(1:2,7)]; Ecc15 = full.expression.path.table[,c(1:2,8)]
E.fae.heatkilled = full.expression.path.table[,c(1:2,9)]; P.rett.heatkilled = full.expression.path.table[,c(1:2,10)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

number.of.cluster=c()
for (g in 1:length(list.of.bacteria.name)){
    sample = data.frame(list.of.bacteria[g]) #Load in the expression path data
    number.of.cluster = c()
    sample[,3] = factor(sample[,3])
    pattern = levels(sample[,3])
    
    for (h in 1:length(pattern)){ #assign genes to clusters
        assign(paste0(list.of.bacteria.name[g],".cluster",h), sample[which(sample[,3] == pattern[g]),])
    }
    cat("Number of clusters for ", list.of.bacteria.name[g], " is ", length(pattern), "\n")
    number.of.cluster[g] = as.integer(length(pattern))
    
    sample.sorted = sample[order(sample[,3]),] #Sort the sample based on expression path and save them
    write.table(sample.sorted, file=paste("edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path_for_",list.of.bacteria.name[g],".txt",sep=""), quote=F, row.names = F, col.names = T )
}

#For a given condition, cluster genes based on patterns (Genes that were EE-EE-EE in ALL conditions were removed)
M.luteus = full.expression.path.table.EEs.removed[,c(1:3)]; E.coli = full.expression.path.table.EEs.removed[,c(1:2,4)]
S.mar.type = full.expression.path.table.EEs.removed[,c(1:2,5)]; E.fae.live = full.expression.path.table.EEs.removed[,c(1:2,6)]
P.rett.live = full.expression.path.table.EEs.removed[,c(1:2,7)]; Ecc15 = full.expression.path.table.EEs.removed[,c(1:2,8)]
E.fae.heatkilled = full.expression.path.table.EEs.removed[,c(1:2,7)]; P.rett.heatkilled = full.expression.path.table.EEs.removed[,c(1:2,8)]
list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), 
                        as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15), as.matrix(E.fae.heatkilled), as.matrix(P.rett.heatkilled))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15", "E.fae.heatkilled", "P.rett.heatkilled")

number.of.cluster.EEs.removed=c()
for (g in 1:length(list.of.bacteria.name)){
    sample = data.frame(list.of.bacteria[g]) #Load in the expression path data
    number.of.cluster.EEs.removed = c()
    sample[,3] = factor(sample[,3])
    pattern = levels(sample[,3])
    
    for (h in 1:length(pattern)){ #assign genes to clusters
        assign(paste0(list.of.bacteria.name[g],".cluster",h), sample[which(sample[,3] == pattern[g]),])
    }
    cat("Number of clusters for ", list.of.bacteria.name[g], " is ", length(pattern), "\n")
    number.of.cluster.EEs.removed[g] = as.integer(length(pattern))
    
    sample.sorted = sample[order(sample[,3]),] #Sort the sample based on expression path and save them
    write.table(sample.sorted, file=paste("edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path_for_",list.of.bacteria.name[g],"_EEs_removed.txt",sep=""), quote=F, row.names = F, col.names = T )
}


#3. Plot the frequency of expression paths for each condition
####################################################
full.expression.path.table.EEs.removed = read.table("edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path_NS_removed.txt", header=T) #1800 x 10
library(plyr)
M.luteus = full.expression.path.table.EEs.removed[,3]; M.luteus.counted = count(M.luteus); M.luteus.counted.sorted = M.luteus.counted[order(M.luteus.counted$x),]
E.coli = full.expression.path.table.EEs.removed[,4]; E.coli.counted = count(E.coli); E.coli.counted.sorted = E.coli.counted[order(E.coli.counted$x),]
S.mar.type = full.expression.path.table.EEs.removed[,5]; S.mar.type.counted = count(S.mar.type); S.mar.type.counted.sorted = S.mar.type.counted[order(S.mar.type.counted$x),]
Ecc15 = full.expression.path.table.EEs.removed[,8]; Ecc15.counted = count(Ecc15); Ecc15.counted.sorted = Ecc15.counted[order(Ecc15.counted$x),]
E.fae.live = full.expression.path.table.EEs.removed[,6]; E.fae.live.counted = count(E.fae.live); E.fae.live.counted.sorted = E.fae.live.counted[order(E.fae.live.counted$x),]
E.fae.heatkilled = full.expression.path.table.EEs.removed[,9]; E.fae.heatkilled.counted = count(E.fae.heatkilled); E.fae.heatkilled.counted.sorted = E.fae.heatkilled.counted[order(E.fae.heatkilled.counted$x),]
P.rett.live = full.expression.path.table.EEs.removed[,7]; P.rett.live.counted = count(P.rett.live); P.rett.live.counted.sorted = P.rett.live.counted[order(P.rett.live.counted$x),]
P.rett.heatkilled = full.expression.path.table.EEs.removed[,10]; P.rett.heatkilled.counted = count(P.rett.heatkilled); P.rett.heatkilled.counted.sorted = P.rett.heatkilled.counted[order(P.rett.heatkilled.counted$x),]

par(mfrow = c(2,4), mar=c(4,4,4,4), oma =c(1,5,0,0))
barplot(M.luteus.counted.sorted$freq, main = "M.luteus", horiz=T, names.arg=M.luteus.counted.sorted$x, las=1)
barplot(E.coli.counted.sorted$freq, main = "E.coli", horiz=T, names.arg=E.coli.counted.sorted$x, las=1)
barplot(S.mar.type.counted.sorted$freq, main = "S.marcescens Type", horiz=T, names.arg=S.mar.type.counted.sorted$x, las=1)
barplot(Ecc15.counted.sorted$freq, main = "Ecc15", horiz=T, names.arg=Ecc15.counted.sorted$x, las=1)
barplot(E.fae.live.counted.sorted$freq, main = "E.faecalis live", horiz=T, names.arg=E.fae.live.counted.sorted$x, las=1)
barplot(E.fae.heatkilled.counted.sorted$freq, main = "E.faecalis heatkilled", horiz=T, names.arg=E.fae.heatkilled.counted.sorted$x, las=1)
barplot(P.rett.live.counted.sorted$freq, main = "P.rettgeri live", horiz=T, names.arg=P.rett.live.counted.sorted$x, las=1)
barplot(P.rett.heatkilled.counted.sorted$freq, main = "P.rettgeri heatkilled", horiz=T, names.arg=P.rett.heatkilled.counted.sorted$x, las=1)

#Remove the genes that had EE-EE-EE in a given condition so that the barplot isn't skewed as much
M.luteus.counted.sorted = M.luteus.counted.sorted[which(M.luteus.counted.sorted$x != 'EE-EE-EE'),]
E.coli.counted.sorted = E.coli.counted.sorted[which(E.coli.counted.sorted$x != 'EE-EE-EE'),]
S.mar.type.counted.sorted = S.mar.type.counted.sorted[which(S.mar.type.counted.sorted$x != 'EE-EE-EE'),]
Ecc15.counted.sorted = Ecc15.counted.sorted[which(Ecc15.counted.sorted$x != 'EE-EE-EE'),]
E.fae.live.counted.sorted = E.fae.live.counted.sorted[which(E.fae.live.counted.sorted$x != 'EE-EE-EE'),]
E.fae.heatkilled.counted.sorted = E.fae.heatkilled.counted.sorted[which(E.fae.heatkilled.counted.sorted$x != 'EE-EE-EE'),]
P.rett.live.counted.sorted = P.rett.live.counted.sorted[which(P.rett.live.counted.sorted$x != 'EE-EE-EE'),]
P.rett.heatkilled.counted.sorted = P.rett.heatkilled.counted.sorted[which(P.rett.heatkilled.counted.sorted$x != 'EE-EE-EE'),]

par(mfrow = c(2,4), mar=c(4,4,4,4), oma =c(1,5,0,0))
barplot(M.luteus.counted.sorted$freq, main = "M.luteus", horiz=T, names.arg=M.luteus.counted.sorted$x, las=1, xlim=c(0,300))
barplot(E.coli.counted.sorted$freq, main = "E.coli", horiz=T, names.arg=E.coli.counted.sorted$x, las=1, xlim=c(0,300))
barplot(S.mar.type.counted.sorted$freq, main = "S.marcescens Type", horiz=T, names.arg=S.mar.type.counted.sorted$x, las=1, xlim=c(0,300))
barplot(Ecc15.counted.sorted$freq, main = "Ecc15", horiz=T, names.arg=Ecc15.counted.sorted$x, las=1, xlim=c(0,300))
barplot(E.fae.live.counted.sorted$freq, main = "E.faecalis live", horiz=T, names.arg=E.fae.live.counted.sorted$x, las=1, xlim=c(0,300))
barplot(E.fae.heatkilled.counted.sorted$freq, main = "E.faecalis heatkilled", horiz=T, names.arg=E.fae.heatkilled.counted.sorted$x, las=1, xlim=c(0,300))
barplot(P.rett.live.counted.sorted$freq, main = "P.rettgeri live", horiz=T, names.arg=P.rett.live.counted.sorted$x, las=1, xlim=c(0,300))
barplot(P.rett.heatkilled.counted.sorted$freq, main = "P.rettgeri heatkilled", horiz=T, names.arg=P.rett.heatkilled.counted.sorted$x, las=1, xlim=c(0,300))


#Create a list of expression paths
expression.list = read.table("edgeR_prev_infected_vs_present_infected_all_genes_with_expression_path.txt", header=T) #11911 x 10

#Pick the longest path
path.list = c()
for (i in c(3:10)){
    path.list = c(path.list, levels(expression.list[,i]))
}
path.list = unique(path.list)
path.list = data.frame(sort(path.list,decreasing = T)) #22 paths
write.table(path.list, file="list_of_expression_paths_from_prev_infected_vs_present_infected.txt", quote=F, row.names = T, col.names = F)






