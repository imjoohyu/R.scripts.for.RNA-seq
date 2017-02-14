#Finding Core Gram-positive and Core Gram-negative genes
#February 11, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/")

data = read.csv("core_genes_all_live_bacteria_upregulated_with_conditions_and_normalized_counts.txt", header=T, sep = "\t")
#data = read.table("core_genes_all_live_bacteria_downregulated_with_conditions_and_normalized_counts.txt", header=T, sep = "\t")
data = data[,1:4]

Gram_positive_exclusive = c(); Gram_negative_exclusive = c()
Core_Gram_positive_exclusive = c(); Core_Gram_negative_exclusive = c()
Shared_Gram_genes = c()

for (i in 1:dim(data)[1]){
    #Format the information on the conditions
    format_information = function(data){
        list_of_conditions = unlist(strsplit(as.character(data[i,4]), ", "))
        list_of_conditions = list_of_conditions[-1]
        
        #strip the time information and only get unique conditions
        for (j in 1:length(list_of_conditions)){
            list_of_conditions[j] = gsub(".12hr","",list_of_conditions[j])
            list_of_conditions[j] = gsub(".36hr","",list_of_conditions[j])
            list_of_conditions[j] = gsub(".5.5d","",list_of_conditions[j])
        }
        return(unique(list_of_conditions))
    }
    list_of_conditions = format_information(data)
    
    #Count the number of conditions
    Gram_positive = c("M.luteus", "E.fae.live", "S.aureus")
    Gram_negative = c("E.coli", "S.mar.type", "P.rett.live", "Ecc15", "P.sneebia", "S.mar.Db11", "P.ento")
    
    count_Gram_occasions = function(Gram_type){
        count = 0
        for (m in 1:length(Gram_type)){
            if (Gram_type[m] %in% list_of_conditions){ #if this condition shows up as DE
                count =  count +1
            }
        }
        return(count)
    }
    Gram_positive_count = count_Gram_occasions(Gram_positive)
    Gram_negative_count = count_Gram_occasions(Gram_negative)
    
    #Make a call based on the answer from above

    
    if (Gram_positive_count > 0 && Gram_negative_count ==0){ #Gram-positive exclusive
        Gram_positive_exclusive = rbind(Gram_positive_exclusive,data[i,1:4])
        
        if (Gram_positive_count == 3){ #Core Gram-positive exclusive
            Core_Gram_positive_exclusive = rbind(Core_Gram_positive_exclusive, data[i,1:4])
        }
    }
    if (Gram_negative_count > 0 && Gram_positive_count ==0){ #Gram-negative exclusive
        Gram_negative_exclusive = rbind(Gram_negative_exclusive,data[i,1:4])
        
        if (Gram_negative_count == 7){ #Core Gram-negative exclusive
            Core_Gram_negative_exclusive = rbind(Core_Gram_negative_exclusive, data[i,1:4])
        }
    }
    else if (Gram_negative_count > 0 & Gram_positive_count >0){
        Shared_Gram_genes = rbind(Shared_Gram_genes, data[i,1:4])
    }
    
}

#Upregulated:
cat("The number of Gram-positive exclusive genes is: ",dim(Gram_positive_exclusive)[1]) #269
cat("The number of Gram-negative exclusive genes is: ",dim(Gram_negative_exclusive)[1]) #464
cat("The number of Core Gram-positive exclusive genes is: ",dim(Core_Gram_positive_exclusive)[1]) #20
cat("The number of Core Gram-negative exclusive genes is: ",dim(Core_Gram_negative_exclusive)[1]) #1
cat("The number of shared genes is: ",dim(Shared_Gram_genes)[1]) #553

#Downregulated:
cat("The number of Gram-positive exclusive genes is: ",dim(Gram_positive_exclusive)[1]) #393
cat("The number of Gram-negative exclusive genes is: ",dim(Gram_negative_exclusive)[1]) #387
cat("The number of Core Gram-positive exclusive genes is: ",dim(Core_Gram_positive_exclusive)[1]) #8
cat("The number of Core Gram-negative exclusive genes is: ",dim(Core_Gram_negative_exclusive)[1]) #NA
cat("The number of shared genes is: ",dim(Shared_Gram_genes)[1]) #510