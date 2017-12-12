#Looking into overal DEG
#March 31st, 2017
#Joo Hyun Im (ji72)

#Identify the number of genes up/downregulated in at least one condition

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")

data = read.table("edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig_direction.txt", header =T) #2589 x 33

data$NumUpDEG = rowSums(data == "Up")
data$NumDownDEG = rowSums(data == "Down")

#Get the number of up/downregulated DEG in at least one bacterial infection (including clean prick and heatkilled)
upDEG_number = dim(data)[1] - table(data$NumUpDEG)[1]; upDEG_number #1374
downDEG_number = dim(data)[1] - table(data$NumDownDEG)[1]; downDEG_number #1386

#Get the number of up/downregulated DEG in at least one bacterial infection (including clean prick and heatkilled)
data_live_only = data[,c(1:2,6:27)]
data_live_only$NumUpDEG = rowSums(data_live_only == "Up")
data_live_only$NumDownDEG = rowSums(data_live_only == "Down")
upDEG_number = dim(data_live_only)[1] - table(data_live_only$NumUpDEG)[1]; upDEG_number #1286
downDEG_number = dim(data_live_only)[1] - table(data_live_only$NumDownDEG)[1]; downDEG_number #1290


#Get the number of DEG that are both up/downregulated
indicator = c()
for (i in 1:dim(data)[1]){
    if (data[i,34] > 0 && data[i,35] > 0){ #if found in both up and down
        indicator = rbind(indicator, c("TRUE")) #both
    }
    else{
        indicator = rbind(indicator, c("FALSE"))
    }
}
data = cbind(data, indicator)
colSums(data == "TRUE") #171


#Get the number of DEG that are DE in 12hr
data_12hr = data[,c(1:2,3,6,9,12,15,18,21,24,25,26,27,28,31)] #all 12hr
data_12hr$NumUpDEG = rowSums(data_12hr == "Up")
data_12hr$NumDownDEG = rowSums(data_12hr == "Down")
upDEG_number_12hr = dim(data_12hr)[1] - table(data_12hr$NumUpDEG)[1]; upDEG_number_12hr #1078
downDEG_number_12hr = dim(data_12hr)[1] - table(data_12hr$NumDownDEG)[1]; downDEG_number_12hr #1092

#Get the number of DEG that are DE in 36hr
data_36hr = data[,c(1:2,4,7,10,13,16,19,22,29,32)]
data_36hr$NumUpDEG = rowSums(data_36hr == "Up")
data_36hr$NumDownDEG = rowSums(data_36hr == "Down")
upDEG_number_36hr = dim(data_36hr)[1] - table(data_36hr$NumUpDEG)[1]; upDEG_number_36hr #493
downDEG_number_36hr = dim(data_36hr)[1] - table(data_36hr$NumDownDEG)[1]; downDEG_number_36hr #314

#Get the number of DEG that are DE in 132hr
data_132hr = data[,c(1:2,5,8,11,14,17,20,23,30,33)]
data_132hr$NumUpDEG = rowSums(data_132hr == "Up")
data_132hr$NumDownDEG = rowSums(data_132hr == "Down")
upDEG_number_132hr = dim(data_132hr)[1] - table(data_132hr$NumUpDEG)[1]; upDEG_number_132hr #501
downDEG_number_132hr = dim(data_132hr)[1] - table(data_132hr$NumDownDEG)[1]; downDEG_number_132hr #260

#The number of upregulated genes at 12h was X times higher than the average number of induced genes at 36h and 132hr
upDEG_number_12hr/((upDEG_number_36hr + upDEG_number_132hr)/2) #2.169x
#The number of downregulated genes at 12h was X times higher than the average number of repressed genes at 36h and 132hr
downDEG_number_12hr/((downDEG_number_36hr + downDEG_number_132hr)/2) #3.804x

#### normalize by only looking at the 6 live conditions
##Get the number of DEG that are DE in 12hr
data_12hr = data[,c(1:2,6,9,12,15,18,21)] #6 live bacteria
data_12hr$NumUpDEG = rowSums(data_12hr == "Up")
data_12hr$NumDownDEG = rowSums(data_12hr == "Down")
data_12hr = data_12hr[which(data_12hr$NumUpDEG != 0 | data_12hr$NumDownDEG != 0),] #1299
upDEG_number_12hr = dim(data_12hr)[1] - table(data_12hr$NumUpDEG)[1]; upDEG_number_12hr #684
downDEG_number_12hr = dim(data_12hr)[1] - table(data_12hr$NumDownDEG)[1]; downDEG_number_12hr #623

#Get the number of DEG that are DE in 36hr
data_36hr = data[,c(1:2,7,10,13,16,19,22)] #6 live bacteria
data_36hr$NumUpDEG = rowSums(data_36hr == "Up")
data_36hr$NumDownDEG = rowSums(data_36hr == "Down")
data_36hr = data_36hr[which(data_36hr$NumUpDEG != 0 | data_36hr$NumDownDEG != 0),] #628
upDEG_number_36hr = dim(data_36hr)[1] - table(data_36hr$NumUpDEG)[1]; upDEG_number_36hr #412
downDEG_number_36hr = dim(data_36hr)[1] - table(data_36hr$NumDownDEG)[1]; downDEG_number_36hr #219

#Get the number of DEG that are DE in 132hr
data_132hr = data[,c(1:2,8,11,14,17,20,23)] #6 live bacteria
data_132hr$NumUpDEG = rowSums(data_132hr == "Up")
data_132hr$NumDownDEG = rowSums(data_132hr == "Down")
data_132hr = data_132hr[which(data_132hr$NumUpDEG != 0 | data_132hr$NumDownDEG != 0),] #668
upDEG_number_132hr = dim(data_132hr)[1] - table(data_132hr$NumUpDEG)[1]; upDEG_number_132hr #449
downDEG_number_132hr = dim(data_132hr)[1] - table(data_132hr$NumDownDEG)[1]; downDEG_number_132hr #219

#The number of upregulated genes at 12h was X times higher than the average number of induced genes at 36h and 132hr
upDEG_number_12hr/((upDEG_number_36hr + upDEG_number_132hr)/2) #1.6x
#The number of downregulated genes at 12h was X times higher than the average number of repressed genes at 36h and 132hr
downDEG_number_12hr/((downDEG_number_36hr + downDEG_number_132hr)/2) #2.8x


#Check the number of DEGs
table(data$log2FC.S.aureus.12hr) #S_aureus_count
table(data$log2FC.P.sneebia.12hr) #P_sneebia_count
1193/187 #6.4x


