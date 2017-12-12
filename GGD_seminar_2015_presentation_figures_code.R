#Figure generator for GGD Wednesday seminar
#March 11th, 2015
#Joo Hyun Im

rm(list=ls(all=TRUE))

#Call the count table of all five samples (aligned by tophat)
total.count <- read.table("/Users/JooHyun/Desktop/HiSeq.totalRNAseq.all/count.data.by.tophat/total.seq.count.table.csv")
total.count <- total.count[1:14869,]

#Figure 1: Some genes are differentially expressed depending on the species and over time.
#Generate two XY plots showing the expression level of each gene at 1) 8hrS vs R, 2) 8hrS vs 24hrS.

#Pull out the data that I want
total.count.subset <- total.count[,2:3]
#total.count.subset <- total.count[,c("count.8hrS","count.24hrS")]
#Normalize the counts relative to each other using DESeq. Change condition names accordingly.
design <- data.frame(row.names = colnames(total.count.subset), condition =c("sneebia","rettgeri"))

#Normalize with DESeq
library("DESeq")
condition = design$condition
cds.all = newCountDataSet(total.count.subset, condition)
#Normalization: estimate the effective library size.
cds.all = estimateSizeFactors(cds.all)
sizeFactors(cds.all)
#Now divide each column of total.count.subset by the size factor
total.count.subset.normalized <- cbind((total.count.subset[1]/sizeFactors(cds.all)[1]),(total.count.subset[2]/sizeFactors(cds.all)[2]))
total.count.subset <-total.count.subset.normalized

#Remove genes whose expression count is under 5
key.mx <-matrix(NA, ncol=2,nrow=dim(total.count)[1])
for (i in 1:dim(total.count)[1]){
    if (total.count.subset[i,1] > 5 && total.count.subset[i,2] >5 ){
        key.mx[i,1] <- total.count.subset[i,1]
        key.mx[i,2] <- total.count.subset[i,2]
    }
}
#get rid of rows with 'NA'
key.mx <- na.omit(key.mx)
x <-log2(key.mx[,1])
y <- log2(key.mx[,2])
#Draw a plot
library(ggplot2)
qplot(x,y,asp=1)+ xlab("log2(counts from 8 hours post-infection)") + ylab("log2(counts from 24 hours post-infection)") + ggtitle("Comparison of gene expression in samples infected with P. sneebia") + theme_classic(base_size = 20) + coord_equal(ratio=1) + geom_abline(intercept=0, slope=1, color="red")
#coord_cartesian(xlim = c(0, 15), ylim =c(0,15)
count.lm = lm(y ~ x)
summary(count.lm)$r.squared

#Figure 2: Some of these genes are significantly differentially expressed (by fold change).
#Generate an MA plot of 8hrS vs 8hrR
total.count.subset <- total.count[,2:3]
#Filter out the ones that had fewer than 5 reads.
filtered <- vector()
for (i in 1:dim(total.count.subset)[1]){
    if (total.count.subset[i,1] > 5 && total.count.subset[i,2] >5 ){
        filtered <- rbind(filtered,total.count.subset[i,1:2]) 
    }
}

total.count.subset <- filtered

#Prep work for DESeq
colnames(total.subset.count) <- c("8hrS","8hrR")
design <- data.frame(row.names = colnames(total.count.subset), condition =c("sneebia","rettgeri"))
condition = design$condition
library("DESeq")
cds = newCountDataSet(total.count.subset, condition)
#Normalization: estimate the effective library size.
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized=TRUE))
cdsB = estimateDispersions(cds, method ="blind", sharingMode="fit-only")
res = nbinomTest(cdsB, "sneebia", "rettgeri")
#plotMA(res)
#To change the ylim to be wider:
plotMA = function(x, ylim,
                  col = ifelse(x$padj>=0.1, "gray32", "red3"),
                  linecol = "#ff000080",
                  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
                  log = "x", cex=0.45, ...)
{
    if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
        stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
    
    x = subset(x, baseMean!=0)
    py = x$log2FoldChange
    if(missing(ylim))
        ylim = c(-2,2) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
    plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
         log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
         cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
    abline(h=0, lwd=4, col=linecol)
}

#Figure 3. DE between sneebia vs rettgeri at 8 hours
#1) Immune genes vs non-immune genes
slices <- c(30,58)
lbels <- c("Immune genes","Non-immune genes")
labels.ingraph=c("34% (30)","66% (58)")
cols=c("orange","steelblue3")
pie(slices, labels = labels.ingraph, col=cols, cex=1.5)
legend("bottom",lbels,cex=1.4,fill=cols, adj=c(0,0.5),x.intersp = 1, y.intersp = 1)

#2) Out of non-immune genes, categories by function
slices <- c(10,4,4,24,16)
lbels <- c("Metabolism","Proteolysis","Development","Others","Unknown")
labels.ingraph=c("17% (10)","7% (4)","7% (4)","41% (24)","28% (16)")
cols=c("lightsteelblue1","lightsteelblue2", "lightsteelblue3","lightsteelblue4", "lightsteelblue")
pie(slices, labels = labels.ingraph, col=cols, cex=1.5)
legend("bottomleft",lbels,cex=1.35, fill=cols, adj=c(0,0.5),x.intersp = 1, y.intersp = 1)

#Figure 4. DESeq interaction plot for all four samples + base-level (0hr)
library("DESeq")
total.count.subset <- total.count[,2:5]
#Filter out the ones that had fewer than 5 reads.
filtered <- vector()
for (i in 1:dim(total.count.subset)[1]){
    if (total.count.subset[i,1] > 5 && total.count.subset[i,2] >5 && total.count.subset[i,3] > 5 && total.count.subset[i,4] >5){
        filtered <- rbind(filtered,total.count.subset[i,1:4])
    }
}

total.count.subset <- filtered

design2 <- data.frame(row.names = colnames(total.count.subset), species = c("S","R","S","R"), time = c("eight","eight","twenty.four","twenty.four"))
cds.all = newCountDataSet(total.count.subset, design2)
#Normalization: estimate the effective library size.
cds.all = estimateSizeFactors(cds.all)
sizeFactors(cds.all)
#Estimate variance without replicates
cdsBlind = estimateDispersions(cds.all, method="blind", sharingMode="fit-only") #if replicates are not available, estimate across conditions.
plotDispEsts(cdsBlind)

#Fit the GLM
fit1.int = fitNbinomGLMs(cdsBlind, count ~ species + time + species*time) #full model
fit0.int = fitNbinomGLMs(cdsBlind, count ~ species + time)
pvalsGLM.int = nbinomGLMTest(fit1.int, fit0.int)
padjGLM.int = p.adjust(pvalsGLM.int, method="BH") 
#Biologically, these are the genes that changes differently between species over time.
res2=cbind(fit1.int, pval=pvalsGLM.int, padj = padjGLM.int)
res.sig2 = res2[which(res$padj <0.1),] #Pull out genes that have a significant species*time interaction effect
res.sig2 = res.sig2[order(res.sig[,4], decreasing=T),] #ranking by the degree of interaction

#Find cases of positive and negative interaction coefficients and then plot an interaction plot 
#including 0hr using normalized counts (counts normalized to size factor)
total.count.subset.for.interaction <- total.count
filtered <- vector()
for (i in 1:dim(total.count.subset.for.interaction)[1]){
    if (total.count.subset.for.interaction[i,1] > 5 && total.count.subset.for.interaction[i,2] >5 && total.count.subset.for.interaction[i,3] > 5 && total.count.subset.for.interaction[i,4] >5 && total.count.subset.for.interaction[i,5] >5){
        filtered <- rbind(filtered,total.count.subset.for.interaction[i,1:5])
    }
}

total.count.subset.for.interaction <- filtered
design.for.interaction <- data.frame(row.names = colnames(total.count.subset.for.interaction), species = c("zero","S","R","S","R"), time = c("zero","eight","eight","twenty.four","twenty.four"))
cds.interaction = newCountDataSet(total.count.subset.for.interaction, design.for.interaction)
#Normalization: estimate the effective library size.
cds.interaction = estimateSizeFactors(cds.interaction)
sizeFactors(cds.interaction)
total.count.subset.normalized <- log10(cbind((total.count.subset.for.interaction[1]/sizeFactors(cds.interaction)[1]),(total.count.subset.for.interaction[2]/sizeFactors(cds.interaction)[2]),(total.count.subset.for.interaction[3]/sizeFactors(cds.interaction)[3]),(total.count.subset.for.interaction[4]/sizeFactors(cds.interaction)[4]),(total.count.subset.for.interaction[5]/sizeFactors(cds.interaction)[5])))
#Cases of positive interaction coefficients
one <- total.count.subset.normalized[c("FBgn0000279"),] #CecC
two <- total.count.subset.normalized[c("FBgn0000276"),] #CecA1
three = total.count.subset.normalized[c("FBgn0038530"),] #AttD
four = total.count.subset.normalized[c("FBgn0004240"),] #Dpt
five = total.count.subset.normalized[c("FBgn0003162"),] #Punch
six = total.count.subset.normalized[c("FBgn0260746"),] #Etc3
seven <- total.count.subset.normalized[c("FBgn0037906"),] #PGRP-LB
eight <- total.count.subset.normalized[c("FBgn0014018"),] #Relish
#Those with negative interaction coefficients
mi1 = total.count.subset.normalized[c("FBgn0087002"),] #Rfabg/Rfabp
mi2 = total.count.subset.normalized[c("FBgn0016075"),] #viking
mi3 = total.count.subset.normalized[c("FBgn0000299"),] #Collagen IV
mi4 = total.count.subset.normalized[c("FBgn0036262"),] #CG6910
mi5 = total.count.subset.normalized[c("FBgn0040099"),] #Lectin-28C

#Put them into a table
#case.studies = rbind(one,two,three,four,seven,eight)
case.studies = rbind(mi4,mi5,six)
#row.names(case.studies) = c("Cecropin C", "Cecropin A1", "Attacin-D", "Diptericin","PGRP-LB","Relish")
row.names(case.studies) = c("CG6910", "Lectin-28C","Ect3")
par(mfrow = c(1,3))
for (i in 1:dim(case.studies)[1]){
    sneebia <- c(case.studies[i,1],case.studies[i,2], case.studies[i,4])
    rettgeri <- c(case.studies[i,1],case.studies[i,3], case.studies[i,5])
    plot(sneebia, type="o", ylim=c(0,5), axes=F,ann=F, col="red")
    title(xlab="Time", ylab="log10(Nor. Read Count)",main=c(rownames(case.studies[i,])), cex.lab=1.4)
    lines(rettgeri, type="o", col="black")
    axis(1, at=1:3, lab=c("0hr","8hrs", "24hrs"), cex.axis=1.3)
    axis(2, las=1)
    box()
}
#plot(sneebia)
#lines(rettgeri, type="o", col="red")
#Just getting legend
par(mfrow = c(1,1))
plot(sneebia)
legend(1,3.7, expression(italic('P.sneebia'), italic('P.rettgeri')), , cex=1.2, col=c("red","black"), pch=21:21, lty=1:1) 
#cf) legend, making it italic:
# legend('bottomright', bty='n', c('C. elegans range', 'Study area'), cex=0.8, fill=c('light gray', 'white'), border=c('black','black')) 

#c("P.sneebia","P.rettgeri")
#Figure 5.Testing the hypothesis:
#Hypothesis 1: I expect the Toll pathway to be activated in response to Gram-negative infection.
toll1 <- total.count.subset.normalized[c("FBgn0010381"),] #drosomycin
toll1.5 <- total.count.subset.normalized[c("FBgn0010385"),] #defensin
toll2 <- total.count.subset.normalized[c("FBgn0035806"),] #PGRP-SD
toll3 <- total.count.subset.normalized[c("FBgn0010441"),] #pelle
toll4 <- total.count.subset.normalized[c("FBgn0030926"),] #persephone
toll5 <- total.count.subset.normalized[c("FBgn0034407"),] #diptericin B

case.studies.toll = rbind(toll1,toll1.5,toll2,toll3,toll4,toll5)
row.names(case.studies.toll) = c("Drosomycin", "Defensin", "PGRP-SD", "Pelle", "Persephone", "Diptericin B")
par(mfrow = c(2,3))
for (i in 1:dim(case.studies.toll)[1]){
    sneebia.toll <- c(case.studies.toll[i,1],case.studies.toll[i,2], case.studies.toll[i,4])
    rettgeri.toll <- c(case.studies.toll[i,1],case.studies.toll[i,3], case.studies.toll[i,5])
    plot(sneebia.toll, type="o", ylim=c(0,5), axes=F,ann=F, col="red")
    title(xlab="Time", ylab="log10(Nor. Read Count)",main=c(rownames(case.studies.toll[i,])))
    lines(rettgeri.toll, type="o", col="black")
    axis(1, at=1:3, lab=c("0hr","8hrs", "24hrs"))
    axis(2, las=1)
    box()
}
#plot(sneebia)
#lines(rettgeri, type="o", col="red")
#legend(1.1,0.9,c("P.sneebia","P.rettgeri"), cex=1.2, col=c("black","red"), pch=21:22, lty=1:2) 

#Hypothesis 2: What are the non-immune genes differentially activated between the two species?
nimm1 <- total.count.subset.normalized[c("FBgn0003162"),] #punch
nimm2 <- total.count.subset.normalized[c("FBgn0000047"),] #act88F
nimm3 <- total.count.subset.normalized[c("FBgn0003162"),] #ect3
case.studies.nimm = rbind(nimm1,nimm2,nimm3)
row.names(case.studies.nimm) = c("Punch", "Act88F", "Ect3")
par(mfrow = c(1,3))
for (i in 1:dim(case.studies.nimm)[1]){
    sneebia.nimm <- c(case.studies.nimm[i,1],case.studies.nimm[i,2], case.studies.nimm[i,4])
    rettgeri.nimm <- c(case.studies.nimm[i,1],case.studies.nimm[i,3], case.studies.nimm[i,5])
    plot(sneebia.nimm, type="o", ylim=c(0,5), axes=F,ann=F, col="red")
    title(xlab="Time", ylab="log10(Nor. Read Count)",main=c(rownames(case.studies.nimm[i,])))
    lines(rettgeri.nimm, type="o", col="black")
    axis(1, at=1:3, lab=c("0hr","8hrs", "24hrs"))
    axis(2, las=1)
    box()
}

#Figure #6 Replot the BL
bl.total <- read.table("/Users/JooHyun/Desktop/BL_total_080814.dat",sep = ",", header=T)
species <- bl.total[,1]
x <- bl.total[,2]
y <- bl.total[,3]
qplot(x,y, aes=(color=bl.total[,1])) + xlab("Time") + ylab("log(CFU/mL)") + ggtitle("Bacterial Load over Time") + theme_classic(base_size = 20) + coord_equal(ratio=1)
