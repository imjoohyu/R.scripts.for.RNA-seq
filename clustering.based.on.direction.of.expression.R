####################################################
#Cluster samples based on the direction of expression
#September 7, 2016
#Joo Hyun Im (ji72)
####################################################

#Create a table of direction of expression for each gene in each condition
#1 = Up, 0 = EE, -1 = Down
#This method reduces the degree of elevation or reduction in expression level because a 2-fold increase is labeled the same was as 4-fold increase, for instance.

#delete any previous input
rm(list=ls(all=TRUE))
#Load the data. Get the samples ready for DE analysis.
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")


########################################################################################################
#Part I. Unchallenged vs Infected dataset
########################################################################################################
#Read in the data file
total.DE.with.name.and.id = read.table("edgeR_unchallenged_vs_infected_all_genes_FC.txt", header=T) #11911 x 64
print("Filtering out the genes that have NA for p-val ..."); total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) 
total.with.names.sig.indicated = read.table("edgeR_unchallenged_vs_infected_all_genes_FC_sig_changed_to_Y_or_N.txt", header=T) #11911 x 64

#1. Process data
####################################################
#Convert 'Up' to 1, 'EE' to 0, 'Down' to -1
print("Converting the FC information to Up, Down, or Same depending on their degree of significance ...")
length.of.table <- as.numeric(dim(total.with.names.sig.indicated)[1]); width.of.table <- as.numeric(dim(total.with.names.sig.indicated)[2])
indicator <- NULL
for (i in 1:length.of.table){ #1, 2, 3, ... 11911
    #cat("The value for k is: ", k, "\n" )
    for (s in seq(3, width.of.table, 2)){ #3, 5, 7, ... 64
        
        if (total.with.names.sig.indicated[i,s+1] == "Y"){ #if significant
            indicator <- isTRUE(total.DE.with.name.and.id[i,s] > 0) #indicator shows that the direction is positive
            #cat("indicator: ", indicator, "\n")
            if (indicator == TRUE) { #If the case is Up-DEG
                total.with.names.sig.indicated[i,s] = 1
                }
            else { #If the caseis Down-DEG
                total.with.names.sig.indicated[i,s] = -1
            }
        }
        else { #if not significant
            total.with.names.sig.indicated[i,s] = 0
        }
        
    }
}
write.table(total.with.names.sig.indicated, file="individual.time.points/edgeR_unchallenged_vs_infected_all_genes_FC_converted_to_numeric_direction.txt", quote = F, row.names = F) #11911 x 64

#Only pull out directions and discard FDR Y/N
total.with.names.sig.indicated.dir.only = total.with.names.sig.indicated[,c(seq(3,63,2))] #all infections including virulent infections that only have 12hr time point. Exclude gene.name and gene.id.


#2. Cluster the conditions (samples) based on direction of expression
####################################################
#Traqnspose the dataframe -- rows: condition, columns: each gene's expression value
total.with.names.sig.indicated.dir.only.t = t(total.with.names.sig.indicated.dir.only)
clusters = hclust(dist(total.with.names.sig.indicated.dir.only.t))
plot(clusters, main="Hierarchical Clustering of DE between unchallenged and infected")








########################################################################################################
#Part II. Previous time point infected  vs Present time point dataset
########################################################################################################
#delete any previous input
rm(list=ls(all=TRUE))
#Load the data. Get the samples ready for DE analysis.
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/")

#Read in the data file
total.DE.with.name.and.id = read.table("edgeR_prev_infected_vs_present_infected_all_genes_FC.txt", header=T) #11911 x 56 -- excluding conditions that only have 12hr time point
print("Filtering out the genes that have NA for p-val ..."); total.DE.with.name.and.id <- na.omit(total.DE.with.name.and.id) 
total.with.names.sig.indicated = read.table("edgeR_prev_infected_vs_present_infected_all_genes_FC_sig_changed_to_Y_or_N.txt", header = T) #11911 x 56

####################################################
#1. Process data
#Convert 'Up' to 1, 'EE' to 0, 'Down' to -1
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
                total.with.names.sig.indicated[i,s] = 1
            }
            else { #If the caseis Down-DEG
                total.with.names.sig.indicated[i,s] = -1
            }
        }
        else { #if not significant
            total.with.names.sig.indicated[i,s] = 0
        }
        
    }
}
write.table(total.with.names.sig.indicated, file="individual.time.points/edgeR_prev_infected_vs_present_infected_all_genes_FC_converted_to_numeric_direction.txt", quote = F, row.names = F) #11911 x 56

#Only pull out directions and discard FDR Y/N
total.with.names.sig.indicated.dir.only = total.with.names.sig.indicated[,c(seq(3,56,2))] #all infections excluding virulent infections that only have 12hr time point. Exclude gene.name and gene.id.


####################################################
#2. Cluster the conditions (samples) based on direction of expression
#Traqnspose the dataframe -- rows: condition, columns: each gene's expression value
total.with.names.sig.indicated.dir.only.t = t(total.with.names.sig.indicated.dir.only)
clusters = hclust(dist(total.with.names.sig.indicated.dir.only.t))
plot(clusters, main="Hierarchical Clustering of DE between n time point and n+1 timpe point")


########################################################################################################
#Part III. Factor analysis?
########################################################################################################
#3. Factor Analysis
fit <- factanal(total.with.names.sig.indicated.dir.only.t, 3, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 

#PCA?
fit = princomp(total.with.names.sig.indicated.dir.only, cor=T)
summary(fit)
loadings(fit)
fit$scores



