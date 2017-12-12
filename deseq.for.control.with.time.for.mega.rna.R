#DESeq Commands for Mega RNA-seq --- testing controls
#Date: September 2, 2015
#ji72

#Open the count table and create the metadata that go with the count table.
setwd("/workdir/ji72/Aligned_Trimmed_Counted_Aug_2015/")
#for (i in 1:105){
#    paste("ID_",i) <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30")
#}

u.zero.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_1_ACATTA_R1.fastq.Aligned.out.sam.count", row.names=1)
u.zero.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_36_TGCTAT_R1.fastq.Aligned.out.sam.count", row.names=1)
u.zero.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_71_GTGTAG_R1.fastq.Aligned.out.sam.count", row.names=1)
u.twelve.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_2_GAACCT_R1.fastq.Aligned.out.sam.count", row.names=1)
u.twelve.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_37_CGCCTG_R1.fastq.Aligned.out.sam.count", row.names=1)
u.twelve.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_72_GGTATA_R1.fastq.Aligned.out.sam.count", row.names=1)
u.thirtysix.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_3_ACAACG_R1.fastq.Aligned.out.sam.count", row.names=1)
u.thirtysix.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_38_CATCTA_R1.fastq.Aligned.out.sam.count", row.names=1)
u.thirtysix.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_73_GATTGT_R1.fastq.Aligned.out.sam.count", row.names=1)
u.fivehalf.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_4_AGTTGA_R1.fastq.Aligned.out.sam.count", row.names=1)
u.fivehalf.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_39_AAGCTC_R1.fastq.Aligned.out.sam.count", row.names=1)
u.fivehalf.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_74_GTGCCA_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.twelve.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_5_AGGCAT_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.twelve.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_40_CAGGAC_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.twelve.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_75_ACCGTG_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.thirtysix.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_6_GAAGTG_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.thirtysix.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_41_CGCAAC_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.thirtysix.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_76_ACGTCT_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.fivehalf.1 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_7_AACAAG_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.fivehalf.2 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5581_2250_23782_HFKCMBGXX_ID_42_CGCGGA_R1.fastq.Aligned.out.sam.count", row.names=1)
cp.fivehalf.3 <- read.table("Aligned_BDGP6.80_Aug_2015_Trimmed_QC30_5582_2250_23783_HG7JWBGXX_ID_77_AGACCA_R1.fastq.Aligned.out.sam.count", row.names=1)

#All together
countTable <- cbind(u.zero.1,u.zero.2,u.zero.3,u.twelve.1,u.twelve.2,u.twelve.3,
                    u.thirtysix.1, u.thirtysix.2, u.thirtysix.3, u.fivehalf.1,u.fivehalf.2, u.fivehalf.3,
                    cp.twelve.1,cp.twelve.2,cp.twelve.3,cp.thirtysix.1,cp.thirtysix.2,cp.thirtysix.3,
                    cp.fivehalf.1,cp.fivehalf.2,cp.fivehalf.3)
colnames(countTable) <- c("u.zero.1","u.zero.2","u.zero.3","u.twelve.1","u.twelve.2","u.twelve.3","u.thirtysix.1","u.thirtysix.2","u.thirtysix.3",
                          "u.fivehalf.1","u.fivehalf.2","u.fivehalf.3","cp.twelve.1","cp.twelve.2","cp.twelve.3","cp.thirtysix.1","cp.thirtysix.2",
                          "cp.thirtysix.3","cp.fivehalf.1","cp.fivehalf.2","cp.fivehalf.3")
design <- data.frame(row.names = colnames(countTable), 
                     treatment = c(rep("Unmolested",12),rep("Clean.prick",9)), 
                     time = c("zero","zero","zero","twelve","twelve","twelve","thirty.six","thirty.six","thirty.six","five.half","five.half","five.half",
                              "twelve","twelve","twelve","thirty.six","thirty.six","thirty.six","five.half","five.half","five.half"),
                     rep = c(rep(1:3,7)))
CountTable <- countTable[(1:(dim(countTable)[1]-5)),]

#By treatment
#Unmolested
countTable.u <- cbind(u.zero.1,u.zero.2,u.zero.3,u.twelve.1,u.twelve.2,u.twelve.3,
                      u.thirtysix.1, u.thirtysix.2, u.thirtysix.3, u.fivehalf.1,u.fivehalf.2, u.fivehalf.3)
colnames(countTable.u) <- c("u.zero.1","u.zero.2","u.zero.3","u.twelve.1","u.twelve.2","u.twelve.3","u.thirtysix.1","u.thirtysix.2","u.thirtysix.3",
                            "u.fivehalf.1","u.fivehalf.2","u.fivehalf.3")
design.u <- data.frame(row.names = colnames(countTable.u), 
                       treatment = c(rep("Unmolested",12)), 
                       time = c("zero","zero","zero","twelve","twelve","twelve","thirty.six","thirty.six","thirty.six","five.half","five.half","five.half"),
                       rep = c(rep(1:3,4)))
CountTable.u <- countTable.u[(1:(dim(countTable.u)[1]-5)),]

#Clean.prick
countTable.cp <- cbind(cp.twelve.1,cp.twelve.2,cp.twelve.3,cp.thirtysix.1,cp.thirtysix.2,cp.thirtysix.3,
                    cp.fivehalf.1,cp.fivehalf.2,cp.fivehalf.3)
colnames(countTable.cp) <- c("cp.twelve.1","cp.twelve.2","cp.twelve.3","cp.thirtysix.1","cp.thirtysix.2",
                          "cp.thirtysix.3","cp.fivehalf.1","cp.fivehalf.2","cp.fivehalf.3")
design.cp <- data.frame(row.names = colnames(countTable.cp), 
                     treatment = c(rep("Clean.prick",9)), 
                     time = c("twelve","twelve","twelve","thirty.six","thirty.six","thirty.six","five.half","five.half","five.half"),
                     rep = c(rep(1:3,3)))
CountTable.cp <- countTable.cp[(1:(dim(countTable.cp)[1]-5)),]


#---
#Filter the hits with less than 5 reads
type.of.table <- CountTable.cp #or CountTable.cp
design <- design.cp #design.cp

filteredCountTable <- vector()
length.of.table <- as.numeric(dim(type.of.table)[1])
width.of.table <- as.numeric(dim(type.of.table)[2])
indicator <- NULL
for (i in 1:length.of.table){
    for (j in 1:width.of.table){
        #cat("ith row: ", i, "jth column: ",j,"value: ",CountTable[i,j], "\n")
        indicator <- isTRUE(type.of.table[i,j] < 5)
        #cat("indicator: ", indicator, "\n")
        if (indicator == TRUE) break
    }
    if (j == width.of.table){
        filteredCountTable <- rbind(filteredCountTable,type.of.table[i,1:width.of.table])   
    }
}

CountTable <- filteredCountTable

library("DESeq")
cds = newCountDataSet(CountTable, design) #design

#Normalization: estimate the effective library size.
cds = estimateSizeFactors(cds)
sizeFactors(cds) #now divide each column of the count table by this size factor
head(counts(cds, normalized=TRUE))

#Diagnostics - PCA
cdsB = estimateDispersions(cds, method ="blind")
vsd = varianceStabilizingTransformation(cdsB)

#For all control samples
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_controls_only_by_treatment.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("treatment"))
dev.off()
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_controls_only_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_controls_only_by_tratment_and_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("treatment","time"))
dev.off()

#For unchallenged controls
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_unmolested_controls_only_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
#For cp controls
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_cleanprick_controls_only_by_time.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time"))
dev.off()
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_PCA_cleanprick_controls_only_by_time_and_rep.pdf", height=5, width=5) 
plotPCA(vsd, intgroup=c("time","rep"))
dev.off()


#Estimate variance using the CR-adjusted profile likelihood dispersion estimates
cds = estimateDispersions(cds, method = "pooled-CR", modelFormula = count ~ time) #if replicates are available. #"pooled-CR"
#Inspect the estimated dispersions
#pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_estimate_of_dispersion_pooled.pdf", height=5, width=5) 
#plotDispEsts(cds)
#dev.off()

#Inference: calling differential expression for controls
fit1.int = fitNbinomGLMs(cds, count ~ time) #full model
fit0.int = fitNbinomGLMs(cds, count ~ 1) #reduced model
pvalsGLM.int = nbinomGLMTest(fit1.int, fit0.int)
padjGLM.int = p.adjust(pvalsGLM.int, method="BH") 
res=cbind(fit1.int, pval=pvalsGLM.int, padj = padjGLM.int)
res.sig = res[which(res$padj <0.1),]
#res.sig = res.sig[order(res.sig[,4], decreasing=T),] #ranking by the degree of interaction
write.table(res.sig, file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/mega_RNA-seq_DESeq_GLM_sig_interactions_with_coefficients_for_unmolested_controls.txt", quote=F)

#Sanity check of inspecting unadjusted p values
pdf(file="/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/DESeq_GLM_unadjusted_pvals_for_unmolested_controls.pdf", height=5, width=5) 
hist(res$pval, breaks=100)
dev.off()


#=========================
#Comparing unchallenged 12hr, 36hr, and 5.5d samples (without unchallenged 0hr) using DESeq
#November 5, 2015
rm(list=ls(all=TRUE))

CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_with_gene_id_and_names.txt", header=T) #11558 x 107
library("DESeq")

#1. Unchallenged 12h, 36h, 5.5d:
CountTable.uc = CountTable[,c(4,5,6,39,40,41,74,75,76)]; rownames(CountTable.uc) = CountTable[,1]
design = data.frame(row.names = colnames(CountTable.uc), condition = c("twelve","thirty.six","five.half","twelve","thirty.six","five.half","twelve","thirty.six","five.half"))

cds = newCountDataSet(CountTable.uc, design) #design
cds = estimateSizeFactors(cds) #Normalization: estimate the effective library size
sizeFactors(cds) #now divide each column of the count table by this size factor
head(counts(cds, normalized=TRUE))

#PCA: This one, unlike the one from edgeR, indicates that one 12hr and one 36hr sample is different. :(
cdsB = estimateDispersions(cds, method ="blind")
vsd = varianceStabilizingTransformation(cdsB); plotPCA(vsd, intgroup=c("condition"))

#DE
cds = estimateDispersions(cds, method = "pooled-CR", modelFormula = count ~ condition) #if replicates are available. #"pooled-CR"
plotDispEsts(cds)
fit1.int = fitNbinomGLMs(cds, count ~ condition) #full model
fit0.int = fitNbinomGLMs(cds, count ~ 1) #reduced model
pvalsGLM.int = nbinomGLMTest(fit1.int, fit0.int)
padjGLM.int = p.adjust(pvalsGLM.int, method="BH") 
res=cbind(fit1.int, pval=pvalsGLM.int, padj = padjGLM.int)
res.sig = res[which(res$padj <0.1),]
#res.sig = res.sig[order(res.sig[,4], decreasing=T),] #ranking by the degree of interaction
write.table(res.sig, file="/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_DESeq_GLM_sig_interactions_with_coefficients_for_unmolested_controls_minus_0hr.txt", quote=F)
#Results: 15 DEGs

#2. Unchallenged 12h, 36h: Are there any differences among 12hr and 36hr?
CountTable.uc = CountTable[,c(4,5,39,40,74,75)]; rownames(CountTable.uc) = CountTable[,1]
design = data.frame(row.names = colnames(CountTable.uc), condition = c("twelve","thirty.six","twelve","thirty.six","twelve","thirty.six"))
cds = newCountDataSet(CountTable.uc, design) #design
cds = estimateSizeFactors(cds) #Normalization: estimate the effective library size
sizeFactors(cds) #now divide each column of the count table by this size factor
head(counts(cds, normalized=TRUE))

#PCA: This one, unlike the one from edgeR, indicates that one 12hr and one 36hr sample is different. :(
cdsB = estimateDispersions(cds, method ="blind")
vsd = varianceStabilizingTransformation(cdsB); plotPCA(vsd, intgroup=c("condition"))

#DE
cds = estimateDispersions(cds, method = "pooled-CR", modelFormula = count ~ condition) #if replicates are available. #"pooled-CR"
plotDispEsts(cds)
fit1.int = fitNbinomGLMs(cds, count ~ condition) #full model
fit0.int = fitNbinomGLMs(cds, count ~ 1) #reduced model
pvalsGLM.int = nbinomGLMTest(fit1.int, fit0.int)
padjGLM.int = p.adjust(pvalsGLM.int, method="BH") 
res=cbind(fit1.int, pval=pvalsGLM.int, padj = padjGLM.int)
res.sig = res[which(res$padj <0.1),]
#Results: 0 DEGs