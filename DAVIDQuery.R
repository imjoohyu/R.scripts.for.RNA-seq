#DAVIDQuery
#October 12th, 2015
#Joo Hyun Im (ji72)

#delete any previous input
rm(list=ls(all=TRUE))
library("DAVIDQuery")

#1. GO analysis on genes that are significantly activated or repressed only in Staph infection.
staph.genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_mega_RNA-seq_Sept_2015/list_of_genes_only_activated_by_S.aureus_but_not_by_others.txt", header=T)
staph.gene.name = staph.genes[,1]; staph.gene.name = as.character(staph.gene.name)
staph.genes.list = c()
for (i in 1:dim(staph.genes)[1]){
    if (i == as.numeric(dim(staph.genes)[1])){
        staph.genes.list = paste(staph.genes.list,staph.gene.name[i], sep="")
    }
    else {
        staph.genes.list = paste(staph.genes.list,staph.gene.name[i],",",sep="")
    }
}
staph.genes.list = as.character(staph.genes.list)

#staph.genes.list = c("FBgn0033366,FBgn0031373,FBgn0034943,FBgn0031860,FBgn0039114")
result = DAVIDQuery(ids = staph.genes.list, type="FLYBASE_GENE_ID", annot="KEGG_PATHWAY", tool="geneReportFull")


result = DAVIDQuery(ids = staph.genes.list, type="UNIPORT_ASSESSION", annot="KEGG_PATHWAY", tool="geneReportFull", verbose =T, writeHTML=T)
    

