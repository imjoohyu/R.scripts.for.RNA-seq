#maSigPro for Multiple Series Time Course Experiment
h#August 1st, 2016
#Joo Hyun Im (ji72)

#delete any previous input
rm(list=ls(all=TRUE))
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/maSigPro/")
library(maSigPro)

#1. Enter data into workspace
data.rnaseq = read.table("Normalized_counts_from_all_genes_with_all_UCs_filtered_by_edgeR.txt", header=T) #11911 genes after normalized and filtered by edgeR
design.rnaseq = read.table("maSigPro_RNA-seq_experimental_design_for_time-course.txt", header=T)
design.masig = make.design.matrix(design.rnaseq, degree= 2) #degree is 2 because there are three time points

#2. Finding significant genes by computing a regression fit for each gene
fit = p.vector(data.rnaseq, design.masig, Q=0.05, MT.adjust ="BH", min.obs=20)
fit$i #returns the number of significant genes
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05) #finding significant differences
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups") #obtaining lists of significant genes
head(sigs$summary); names(sigs$sig.genes); names(sigs$sig.genes$ColdvsControl)

#3. Visualize the results
#i) Create a venn diagram of the DE genes for some treatment groups
#For instance, look for genes that are DE in both M.luteus and E.coli
suma2Venn(sigs$summary[, c(3:4)]) #graphical display -- venn diagram of the DE genes for all the treatment groups
#ii) Visualize the significant genes of a particular comparison by clustering
#For instance, comparison between P.rettgeri vs control
see.genes(sigs$sig.genes$M.luteusvscontrol, main = "M.luteus vs control", show.fit = T, dis =design.masig$dis, cluster.method="kmeans", cluster.data = 1, k = 9)

#Problem: commented on 8/1/2016
#When plotting see.genes() command, I get graphs that "evaluate the consistency of the clusters" correctly, but the next set of graphs that are supposed to show the differences between groups do not seem to work. I only get such graph for cluster 1 while there're 9 clusters.


