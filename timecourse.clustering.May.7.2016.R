#Clustering analysis of RNA-seq data for kinetics/time-course
#May 7, 2016
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE))

#################################################################################################################
##############################################
#Option 1: Use MBCluster.seq (Model based clustering)
##############################################
#################################################################################################################
#5/9/2016:
#According to the Statstical Analysis of Next Generation Sequencing book, Ch10
"the count data that arise from RNA-seq experiments are far from normally distributed for most
genes. Another problem is that many genes may have low counts or zero counts
in some treatment groups, and this introduces problems in the log transformation.
Alternatively, model-based approaches using Poisson or negative binomial models
can handle this problem easily. Recently, there have been methods proposed for
clustering count data based on Poisson and negative binomial models [37, 43], and
these methods will be described in this subsection."

#According to the Statstical Analysis of Next Generation Sequencing book,
#For option 1: Model-based clustering/K-mean clustering
#First, "apply the MB-EM algorithm to the count data based on mixture of negative binomial models".
#Running this algorithm with k = 2, 3, 4,....50 will lead to AIC reaching its minimum at some number of k.
#This k number will be the number of clusters that we will ask the program to make.

#On the CBSU machine
#1. I will use the package called MBCluster.Seq ###Run this on CBSU
rm(list=ls(all=TRUE))
setwd("/workdir/ji72/")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
CountTable = CountTable[,c(2:35, 37:70, 72:105)] #excluding unchallenged 0hr samples, 17558 x 102
colnames(CountTable) = c("UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
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

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

design = read.table("mega_RNA-seq_experimental_design_for_all_samples.txt", header=T)
unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d) #M.luteus
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d) #E.coli
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d) #S.mar.type
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d) #E.fae.live
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d) #P.rett.live
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d) #Ecc15

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

library(MBCluster.Seq)
for (m in 1:length(list.of.bacteria.name)){
    condition.count.table = list.of.bacteria[[m]] #Determine the number of clusters
    
    if (m == 2){
        mydata = RNASeq.Data(condition.count.table,  Normalize=NULL, Treatment = c("UC","UC","UC","Twelve","Twelve","Twelve","Thirty-six","Thirty-six","Thirty-six","Five.half","Five.half","Five.half"), GeneID=row.names(condition.count.table))
    }
    else {
        mydata = RNASeq.Data(condition.count.table,  Normalize=NULL, Treatment = c("UC","UC","UC","Twelve","Twelve","Twelve","Thirty-six","Thirty-six","Five.half","Five.half","Five.half"), GeneID=row.names(condition.count.table))
    }
    c0 = KmeansPlus.RNASeq(mydata, nK=50,  model ="nbinom")$centers  ## choose 50 cluster centers to initialize the clustering (Takes about 30 minutes)
    write.table(c0, paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[m],"_c0_k_50_clustering.txt", sep=""), quote=F, col.names = T, row.names = F)
    
    #Calculate the log-likelihood of the model and the associated AIC value.
    lglk=c(); aic = c() 
    for (n in 2:50){
        cls=Cluster.RNASeq(data=mydata, model="nbinom", centers=c0[1:n,], method="EM")$cluster ## use EM algorithm to cluster genes #The output 'cls' stores the cluster IDs for each gene.
        lglk.new=lglk.cluster(mydata,model="nbinom",cluster=cls) #calculates the log-likelihood given the clustering results from cls. The output can be used to calculate AIC values to select the number of clusters,k.
        lglk = data.frame(rbind(lglk, lglk.new))
        aic.new = as.numeric((2*4) - 2*(lglk.new)) #The smaller AIC, the better. n = 4? (four time points)
        aic = data.frame(rbind(aic, aic.new))
    }
    write.table(aic, file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[m], "_c0_k_50_clustering_AIC.txt", sep=""), quote=F, row.names = F, col.names = F)
    
    #Plot lglk.aic where the x-axis is number of clusters (k) and the y-axis is AIC. Figure out where the plateau is.
    aic.with.label = cbind(c(2:50),aic); colnames(aic.with.label)  = c("nK","AIC")
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[m],"_clustering_k_50_AIC_plot.pdf", sep=""))
    plot(aic.with.label$AIC~aic.with.label$nK); dev.off()
    
    #The nK that has the minimum AIC is the K. Set that number as the number of clusters.
    minimum = order(aic.with.label[,2])[1]
    cls=Cluster.RNASeq(data=mydata, model="nbinom", centers=c0[1:minimum,], method="EM")$cluster
    tr=Hybrid.Tree(data=mydata,cluster=cls, model="nbinom") ## bulild a tree structure for the resulting k clusters
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[m], "_clustering.pdf", sep=""))
    plotHybrid.Tree(merge=tr, cluster=cls, logFC=mydata$logFC, tree.title=NULL, colorful = TRUE) ## plot the tree structure
    dev.off()
    condition.count.table.with.clustering = cbind(condition.count.table, cls)
    write.table(condition.count.table.with.clustering, file = paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[m], "_clustering_values.txt", sep=""), quote =F, row.names=T, col.names=T)
} #Letting the program decide how many clusters to have

#2. Draw the representation plots and heatmap ****on a local machine****
rm(list=ls(all=TRUE))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
library(ggplot2); library(grid); library(gridExtra)

# # Multiple plot function (from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
# # Not used as of 5/16/2016
# #
# # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# # - cols:   Number of columns in layout
# # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# #
# # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# # then plot 1 will go in the upper left, 2 will go in the upper right, and
# # 3 will go all the way across the bottom.
# #
# #multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#     library(grid)
#     
#     # Make a list from the ... arguments and plotlist
#     plots <- c(list(...), plotlist)
#     
#     numPlots = length(plots)
#     
#     # If layout is NULL, then use 'cols' to determine layout
#     if (is.null(layout)) {
#         # Make the panel
#         # ncol: Number of columns of plots
#         # nrow: Number of rows needed, calculated from # of cols
#         layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                          ncol = cols, nrow = ceiling(numPlots/cols))
#     }
#     
#     if (numPlots==1) {
#         print(plots[[1]])
#         
#     } else {
#         # Set up the page
#         grid.newpage()
#         pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#         
#         # Make each plot, in the correct location
#         for (i in 1:numPlots) {
#             # Get the i,j matrix positions of the regions that contain this subplot
#             matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#             
#             print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                             layout.pos.col = matchidx$col))
#         }
#     }
# }

#Plot clustering results from model-based methods using the negative binomial model and K=50.
for (s in 1:length(list.of.bacteria.name)){
    bacteria.w.clustering.value = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[s], "_clustering_values.txt", sep=""), header=T)
    
    if (s == 2){ #for E.coli missing one of the 36hr replicates
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:5],1,mean),apply(bacteria.w.clustering.value[,6:8],1,mean),apply(bacteria.w.clustering.value[,9:11],1,mean))
        mean.bacteria.log10 = log10(mean.bacteria + 1)
        mean.bacteria.log10.w.value = as.data.frame(cbind(mean.bacteria.log10, bacteria.w.clustering.value[,12]))
    }
    else {
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:6],1,mean),apply(bacteria.w.clustering.value[,7:9],1,mean),apply(bacteria.w.clustering.value[,10:12],1,mean))
        mean.bacteria.log10 = log10(mean.bacteria + 1)
        mean.bacteria.log10.w.value = as.data.frame(cbind(mean.bacteria.log10, bacteria.w.clustering.value[,13]))
    }
    colnames(mean.bacteria.log10.w.value) = c("0hr","12hr","36hr","5.5d","cluster_number")
    zero = data.frame(mean.bacteria.log10.w.value[,1]); twelve = data.frame(mean.bacteria.log10.w.value[,2]); thirty.six = data.frame(mean.bacteria.log10.w.value[,3]); five.half = data.frame(mean.bacteria.log10.w.value[,4])
    zero.tag = cbind(zero, rep("0hr",dim(zero)[1]), mean.bacteria.log10.w.value[,5]); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]), mean.bacteria.log10.w.value[,5])
    ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1]), mean.bacteria.log10.w.value[,5]); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]), mean.bacteria.log10.w.value[,5])
    colnames(zero.tag) = c("log10Exp", "Time","Cluster");colnames(twelve.tag) = c("log10Exp", "Time","Cluster"); colnames(ts.tag) = c("log10Exp", "Time","Cluster"); colnames(fs.tag) = c("log10Exp", "Time","Cluster")
    total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.bacteria.log10.w.value),4)); total = cbind(total,gene.names)
    total$Time <- factor(total$Time); total$Cluster = factor(total$Cluster); total$gene.names = factor(total$gene.names)
    
    #Plot 1: First 15 clusters
    c1 = total[which(total$Cluster == 1),]; c2 = total[which(total$Cluster == 2),];c3 = total[which(total$Cluster == 3),]
    c4 = total[which(total$Cluster == 4),]; c5 = total[which(total$Cluster == 5),];c6 = total[which(total$Cluster == 6),]
    c7 = total[which(total$Cluster == 7),]; c8 = total[which(total$Cluster == 8),];c9 = total[which(total$Cluster == 9),]
    c10 = total[which(total$Cluster == 10),]; c11 = total[which(total$Cluster == 11),];c12 = total[which(total$Cluster == 12),]
    c13 = total[which(total$Cluster == 13),]; c14 = total[which(total$Cluster == 14),];c15 = total[which(total$Cluster == 15),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_first_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()
    
    #Plot 2: Biggest 15 clusters (This isn't very useful because there are too many entries)
    order.of.c.size = order(summary(total$Cluster),decreasing =T)
    c1 = total[which(total$Cluster == order.of.c.size[1]),]; c2 = total[which(total$Cluster == order.of.c.size[2]),];c3 = total[which(total$Cluster == order.of.c.size[3]),]
    c4 = total[which(total$Cluster == order.of.c.size[4]),]; c5 = total[which(total$Cluster == order.of.c.size[5]),];c6 = total[which(total$Cluster == order.of.c.size[6]),]
    c7 = total[which(total$Cluster == order.of.c.size[7]),]; c8 = total[which(total$Cluster == order.of.c.size[8]),];c9 = total[which(total$Cluster == order.of.c.size[9]),]
    c10 = total[which(total$Cluster == order.of.c.size[10]),]; c11 = total[which(total$Cluster == order.of.c.size[11]),]; c12 = total[which(total$Cluster == order.of.c.size[12]),]
    c13 = total[which(total$Cluster == order.of.c.size[13]),]; c14 = total[which(total$Cluster == order.of.c.size[14]),]; c15 = total[which(total$Cluster == order.of.c.size[15]),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()

    #     plot.list = list()
#     for (p in 1:15){
#         c = total[which(total$Cluster == order.of.c.size[p]),]
#         plot.ind = ggplot(data =c, aes(x=c$Time, y = c$log10Exp, group = c$gene.names)) + geom_line() + theme_classic()
#         plot.list[[p]] = plot.ind
#     }
#     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
#     multiplot(plot.list[[1]], cols=5)
#     dev.off()
    
}

#leftover code (5/15/2016)
# plot.list = list()
# for (p in 1:15){
#     mean.M.luteus.subset = mean.M.luteus.log10.w.value[which(mean.M.luteus.log10.w.value$cluster_number == p),]
#     mean.M.luteus.subset = mean.M.luteus.subset[,1:4]
#     zero = data.frame(mean.M.luteus.subset[,1]); twelve = data.frame(mean.M.luteus.subset[,2]); thirty.six = data.frame(mean.M.luteus.subset[,3]); five.half = data.frame(mean.M.luteus.subset[,4])
#     zero.tag = cbind(zero, rep("0hr",dim(zero)[1])); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]))
#     ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1])); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]))
#     colnames(zero.tag) = c("v1", "v2");colnames(twelve.tag) = c("v1", "v2");colnames(ts.tag) = c("v1", "v2");colnames(fs.tag) = c("v1", "v2");
#     total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.M.luteus.subset),4))
#     total = cbind(total,gene.names)
#     total$v2 <- factor(total$v2); total$gene.names <- factor(total$gene.names)
#     plot(total$v2, total$v1)
#     #plot = print(ggplot(data = total, aes(x=total$v2, y=total$v1, group=total$gene.names)) + geom_line())
# }
# mean.M.luteus.subset = mean.M.luteus.log10.w.value[which(mean.M.luteus.log10.w.value$cluster_number == p),]
# mean.M.luteus.subset = mean.M.luteus.subset[,1:4]
# plot = ggplot(data = total, aes(x=total$Time, y=total$log10Exp, group=total$gene.names)) + geom_line()  #group=total$gene.names, 
# plot + facet_wrap(~Cluster, nrow=3) + theme_classic()



#################################################################################################################
##############################################
#Option 2: Use regular cluster package for kmeans and hierarchical clustering
##############################################
#################################################################################################################

#1. Use DESeq2 to transform data #Run this on CBSU
library("DESeq2"); setwd("/workdir/ji72")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt"); dim(CountTable) #17558 x 105

#Set the tags
count.table <- cbind(CountTable[,c(2:35,37:70,72:105)]) #everything excluding 0hr unmolested
list.three = c(rep(c(  rep("UC",each=3), rep("Clean.prick",each=3), rep("M.luteus",each=3), rep("E.coli",each=3), 
                       rep("S.mar.type",each=3), rep("E.faecalis",each=3), rep("P.rettgeri",each=3), rep("Ecc15",each=3), 
                       "S.aureus", "P.sneebia", "S.mar.DB11","P.entomophila", rep("E.faecalis.heatkilled",each=3), 
                       rep("P.rettgeri.heatkilled",each=3)),3))
time.three = c(rep( c(c("zero", "zero", "zero"), rep(c("twelve","thirty.six","five.half"), 7),rep("twelve", 4),rep(c("twelve","thirty.six","five.half"),2)),3) )

design = data.frame(row.names = colnames(count.table), treatment = list.three, time = time.three, rep = c(rep(1:3, each = 34)) )
group = factor(paste0(design$treatment, ".", design$time))

cds = DESeqDataSetFromMatrix(countData = count.table, colData = design, design = ~ treatment)
cds = DESeq(cds) #estimating size factors, estimating dispersions, gene-wise dispersion estimates, mean-dispersion relationship, final dispersion estimates
vsd <- varianceStabilizingTransformation(cds, blind=FALSE) #17558 x 102

head(assay(vsd), 3)
vsd.table = assay(vsd)
write.table(vsd.table, file="/home/ji72/RNAseq/mega_RNA-seq_count_of_all_samples_Apr_2016_var_sta_transformed_May_2016.txt", quote=F, col.names = T, row.names = F)

#Comment: VST produces transformed data on the log2 scale which has been normalized with respect to library size.
#The point is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm 
#of count data when the mean is low. (As mean increases, variance also increases in RNA-seq data).
#It uses the experiment-wide trend of variance over mean, in order to transform the data to remove the experiment-wide trend. 
#Note that we do not require or desire that all the genes have exactly the same variance after transformation 
#From DESeq2 manual: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

#2. Run cluster R package on the CBSU machine - VST
rm(list=ls(all=TRUE))
vsd.table = read.table("/workdir/ji72/mega_RNA-seq_count_of_all_samples_Apr_2016_var_sta_transformed_May_2016.txt", header=T)
gene.name = read.table("/workdir/ji72/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
row.names(vsd.table) = row.names(gene.name)

colnames(vsd.table) = c("UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
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

unchallenged.0hr = vsd.table[,c(1:3,35:37,69:71)]; clean.prick.12hr = vsd.table[,c(4,38,72)]; clean.prick.36hr = vsd.table[,c(5,39,73)];clean.prick.5.5d = vsd.table[,c(6,40,74)]; M.luteus.12hr = vsd.table[,c(7,41,75)]; M.luteus.36hr = vsd.table[,c(8,42,76)]; M.luteus.5.5d = vsd.table[,c(9,43,77)]; E.coli.12hr = vsd.table[,c(10,44,78)]; E.coli.36hr = vsd.table[,c(11,79)]; 
E.coli.5.5d = vsd.table[,c(12,46,80)]; S.mar.type.12hr = vsd.table[,c(13,47,81)]; S.mar.type.36hr = vsd.table[,c(14,48,82)]; S.mar.type.5.5d = vsd.table[,c(15,49,83)]; E.fae.live.12hr = vsd.table[,c(16,50,84)]; E.fae.live.36hr = vsd.table[,c(17,51,85)]; E.fae.live.5.5d = vsd.table[,c(18,52,86)]; P.rett.live.12hr = vsd.table[,c(19,53,87)];
P.rett.live.36hr = vsd.table[,c(20,54,88)]; P.rett.live.5.5d = vsd.table[,c(21,55,89)]; Ecc15.12hr = vsd.table[,c(22,56,90)]; Ecc15.36hr = vsd.table[,c(23,57,91)]; Ecc15.5.5d = vsd.table[,c(24,58,92)]; S.aureus.12hr = vsd.table[,c(25,59,93)]; P.sneebia.12hr = vsd.table[,c(26,60,94)]; S.mar.Db11.12hr = vsd.table[,c(27,61,95)];
P.ento.12hr = vsd.table[,c(28,62,96)]; E.fae.heatkilled.12hr = vsd.table[,c(29,63,97)]; E.fae.heatkilled.36hr = vsd.table[,c(30,64,98)]; E.fae.heatkilled.5.5d = vsd.table[,c(31,65,99)]; P.rett.heatkilled.12hr = vsd.table[,c(32,66,100)]; P.rett.heatkilled.36hr = vsd.table[,c(33,67,101)]; P.rett.heatkilled.5.5d = vsd.table[,c(34,68,102)]

unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d) #M.luteus
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d) #E.coli
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d) #S.mar.type
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d) #E.fae.live
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d) #P.rett.live
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d) #Ecc15

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

#Method from: http://www.statmethods.net/advstats/cluster.html
library(cluster) 
for (t in 1:length(list.of.bacteria.name)){
    # I. Partitioning
    condition.count.table = list.of.bacteria[[t]] #Determine the number of clusters
    wss <- (nrow(condition.count.table)-1)*sum(apply(condition.count.table,2,var)) 
    for (i in 2:15) {wss[i] <- sum(kmeans(condition.count.table,centers=i)$withinss) }
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_WSS_based_on_VST.pdf", sep =""))
    plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") #the number plateaus at about 10
    dev.off()
    
    # II. K-Means Cluster Analysis
    fit <- kmeans(condition.count.table, 10) # 10 cluster solution
    aggregate(condition.count.table, by=list(fit$cluster),FUN=mean) # get cluster means 
    condition.count.table <- cbind(condition.count.table, fit$cluster) # append cluster assignment
    colnames(condition.count.table)[ncol(condition.count.table)] <- "cluster_info"
    # vary parameters for most readable graph
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_cluster_based_on_VST.pdf", sep=""))
    clusplot(condition.count.table, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
    dev.off()
    
    # III. Ward Hierarchical Clustering
    d <- dist(condition.count.table, method = "euclidean") # distance matrix
    fit <- hclust(d, method="ward")
    groups <- cutree(fit, k=10) # cut tree into 10 clusters
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_hierarchical_cluster_based_on_VST.pdf", sep=""))
    plot(fit) # display dendogram
    rect.hclust(fit, k=10, border="red")# draw dendogram with red borders around the 10 clusters
    dev.off()
    
    write.table(condition.count.table, file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[t], "_kmeans_nK10_clustering_info_based_on_VST.txt", sep=""), quote=F, row.names = T, col.names = T)
}

#3. Draw the representation plots on a local machine
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
library(ggplot2); library(grid); library(gridExtra)
for (s in 1:length(list.of.bacteria.name)){
    bacteria.w.clustering.value = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.kmeans.based.on.VST/", list.of.bacteria.name[s], "_kmeans_nK10_clustering_info_based_on_VST.txt", sep=""), header=T)
    
    if (s == 2){ #for E.coli missing one of the 36hr replicates
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:5],1,mean),apply(bacteria.w.clustering.value[,6:8],1,mean),apply(bacteria.w.clustering.value[,9:11],1,mean),bacteria.w.clustering.value[,12])
    }
    else {
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:6],1,mean),apply(bacteria.w.clustering.value[,7:9],1,mean),apply(bacteria.w.clustering.value[,10:12],1,mean),bacteria.w.clustering.value[,13])
    }
    colnames(mean.bacteria) = c("0hr","12hr","36hr","5.5d","cluster_number")
    zero = data.frame(mean.bacteria[,1]); twelve = data.frame(mean.bacteria[,2]); thirty.six = data.frame(mean.bacteria[,3]); five.half = data.frame(mean.bacteria[,4])
    zero.tag = cbind(zero, rep("0hr",dim(zero)[1]), mean.bacteria[,5]); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]), mean.bacteria[,5])
    ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1]), mean.bacteria[,5]); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]), mean.bacteria[,5])
    colnames(zero.tag) = c("log10Exp", "Time","Cluster");colnames(twelve.tag) = c("log10Exp", "Time","Cluster"); colnames(ts.tag) = c("log10Exp", "Time","Cluster"); colnames(fs.tag) = c("log10Exp", "Time","Cluster")
    total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.bacteria),4)); total = cbind(total,gene.names)
    total$Time <- factor(total$Time); total$Cluster = factor(total$Cluster); total$gene.names = factor(total$gene.names)
    
    c1 = total[which(total$Cluster == 1),]; c2 = total[which(total$Cluster == 2),];c3 = total[which(total$Cluster == 3),]
    c4 = total[which(total$Cluster == 4),]; c5 = total[which(total$Cluster == 5),];c6 = total[which(total$Cluster == 6),]
    c7 = total[which(total$Cluster == 7),]; c8 = total[which(total$Cluster == 8),];c9 = total[which(total$Cluster == 9),]
    c10 = total[which(total$Cluster == 10),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.kmeans.based.on.VST/",list.of.bacteria.name[s],"_kmeans_nK10_clustering_info_based_on_VST_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10, ncol=5)
    dev.off()
}


#################################################################################################################
##############################################
#Option 3: Only use the genes that were DE in all 3 time points (12hr, 36hr, and 5.5d samples).
#DE between each of these time points and the 0hr unchallenged samples. This way may show a stronger impact of infection
#After pulling out the DE genes, draw the clustering heatmap using the [EBSeq-HMM] method.
##############################################
#################################################################################################################

#5/16/2016
de.list = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)

#(i) Separate the table into 10 different conditions and pull out only the genes that are DE in all time points
de.list = de.list[,c(1,2,9:52)]
M.luteus = de.list[,c(1:2,3:8)]; E.coli = de.list[,c(1:2,9:14)]; S.mar.type = de.list[,c(1:2,15:20)]; E.fae.live = de.list[,c(1:2,21:26)]; P.rett.live =de.list[,c(1:2,27:32)]; Ecc15 = de.list[,c(1:2,33:38)]

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

for (i in 1:length(list.of.bacteria.name)){
    condition.count.table = list.of.bacteria[[i]]
    sig.DE.only = condition.count.table[!grepl("N",condition.count.table[,4]),]
    sig.DE.only = sig.DE.only[!grepl("N",sig.DE.only[,6]),]
    sig.DE.only = sig.DE.only[!grepl("N",sig.DE.only[,8]),]
    write.table(sig.DE.only, file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/",list.of.bacteria.name[i],"_sig_DE_in_all_time_points.txt",sep=""), quote=F, row.names = F, col.names = T)
}

#(ii) Here's the list of DE genes that are significant in all time points in each infection condition
M.luteus.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/M.luteus_sig_DE_in_all_time_points.txt", header=T)
M.luteus.gene.only = M.luteus.gene.only[,1]
E.coli.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.coli_sig_DE_in_all_time_points.txt", header=T)
E.coli.gene.only = E.coli.gene.only[,1]
S.mar.type.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/S.mar.type_sig_DE_in_all_time_points.txt", header=T)
S.mar.type.gene.only = S.mar.type.gene.only[,1]
E.fae.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/E.fae.live_sig_DE_in_all_time_points.txt", header=T)
E.fae.live.gene.only = E.fae.live.gene.only[,1]
P.rett.live.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/P.rett.live_sig_DE_in_all_time_points.txt", header=T)
P.rett.live.gene.only = P.rett.live.gene.only[,1]
Ecc15.gene.only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.only.DE.genes/Ecc15_sig_DE_in_all_time_points.txt", header=T)
Ecc15.gene.only = Ecc15.gene.only[,1]

#(iii) Create a heatmap in a local machine
CountTable = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
CountTable = CountTable[,c(2:35, 37:70, 72:105)] #excluding unchallenged 0hr samples, 17558 x 102
colnames(CountTable) = c("UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
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

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d); M.luteus = M.luteus[row.names(M.luteus) %in% M.luteus.gene.only,]
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d); E.coli = E.coli[row.names(E.coli) %in% E.coli.gene.only,]
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d); S.mar.type = S.mar.type[row.names(S.mar.type) %in% S.mar.type.gene.only,]
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d); E.fae.live = E.fae.live[row.names(E.fae.live) %in% E.fae.live.gene.only,]
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d); P.rett.live = P.rett.live[row.names(P.rett.live) %in% P.rett.live.gene.only,]
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d); Ecc15 = Ecc15[row.names(Ecc15) %in% Ecc15.gene.only,]

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

for (j in 1:length(list.of.bacteria.name)){
    cat("EBSeq-HMM - We are working on the following comparison: 0hr-12hr-36hr-5.5d of ", as.character(list.of.bacteria.name[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.bacteria[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 2){ #getting rid of Rep2 of E.coli 36hr
        CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d")} #removing 36hr rep2 sample 
    else {CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d") }    
    Conditions <- factor(CondVector, levels=c("t.0hr","t.12hr","t.36hr","t.5.5d"))
    
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=count.table, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
    EBSeqHMMGeneOut$MgAllMAPChar #classifies the pattern of gene expression over time per gene
    AllPaths = GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=TRUE) #Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
    gene.list.for.a.given.path = GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, cutoff=0.8, OnlyDynamic=TRUE,Paths=NULL) #Obtain confident gene calls for classifying genes into expression paths, PP>=0.5 is most likely path
    number.of.genes.in.each.category = gene.list.for.a.given.path$NumEach; number.of.genes.in.each.category #already sorted
    genes.in.each.category = gene.list.for.a.given.path$EachPathNames; genes.in.each.category
    de.list = GetDECalls(EBSeqHMMGeneOut,FDR=.05) #Another way to obtain DE gene/isoform list at a certain FDR (but not PPs)
    
    #For a specific gene, draw a plot of expression
    #PlotExp(count.table, Conditions, "Gene_1")
    
    write.table(de.list, file=paste("EBSeqHMM_DE_results_for_",list.of.bacteria.name[j],".txt", sep=""), quote=F)
}


#Heatmap of (log10 + 1) normalization
library("ggplot2")
M.luteus.nor = as.matrix(log10(M.luteus)+1); is.na(M.luteus.nor) <- sapply(M.luteus.nor, is.infinite)
heatmap.2(M.luteus.nor, distfun = dist, hclustfun = hclust, dendrogram=c("row"), col="cm.colors")
heatmap(M.luteus.nor, distfun = dist, hclustfun = hclust, dendrogram=c("row"), col="cm.colors")

#Or heatmap based on fold change
yala = as.matrix(M.luteus.gene.only[,c(3,5,7)]); row.names(yala) = M.luteus.gene.only[,1]
heatmap.2(yala, distfun = dist, hclustfun = hclust, dendrogram=c("row"), col="cm.colors")
yala = as.matrix(E.coli.gene.only[,c(3,5,7)]); row.names(yala) = E.coli.gene.only[,1]
heatmap.2(yala, distfun = dist, hclustfun = hclust, dendrogram=c("row"), col="cm.colors")


#################################################################################################################
##############################################
#Option 4: Do the clustering using MBCluster.seq with nK=27 (3 x 3 x 3)
##############################################
#################################################################################################################

#1. I will use the package called MBCluster.Seq ###Run this on CBSU
rm(list=ls(all=TRUE))
setwd("/workdir/ji72/")
CountTable = read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
CountTable = CountTable[,c(2:35, 37:70, 72:105)] #excluding unchallenged 0hr samples, 17558 x 102
colnames(CountTable) = c("UC.Rep1","UC.Rep2","UC.Rep3","clean.prick.12hr.Rep1","clean.prick.36hr.Rep1","clean.prick.5.5d.Rep1",
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

unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]

unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d) #M.luteus
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d) #E.coli
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d) #S.mar.type
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d) #E.fae.live
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d) #P.rett.live
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d) #Ecc15

list.of.bacteria = list(as.matrix(M.luteus), as.matrix(E.coli), as.matrix(S.mar.type), as.matrix(E.fae.live), as.matrix(P.rett.live), as.matrix(Ecc15))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")

library(MBCluster.Seq)
for (m in 1:length(list.of.bacteria.name)){
    condition.count.table = list.of.bacteria[[m]] #Determine the number of clusters
    
    if (m == 2){
        mydata = RNASeq.Data(condition.count.table,  Normalize=NULL, Treatment = c("UC","UC","UC","Twelve","Twelve","Twelve","Thirty-six","Thirty-six","Thirty-six","Five.half","Five.half","Five.half"), GeneID=row.names(condition.count.table))
    }
    else {
        mydata = RNASeq.Data(condition.count.table,  Normalize=NULL, Treatment = c("UC","UC","UC","Twelve","Twelve","Twelve","Thirty-six","Thirty-six","Five.half","Five.half","Five.half"), GeneID=row.names(condition.count.table))
    }
    c0 = KmeansPlus.RNASeq(mydata, nK=27,  model ="nbinom")$centers  ## choose 27 cluster centers to initialize the clustering
    write.table(c0, paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/", list.of.bacteria.name[m],"_c0_k_27_clustering.txt", sep=""), quote=F, col.names = T, row.names = F)
    cls=Cluster.RNASeq(data=mydata, model="nbinom", centers=c0[1:27,], method="EM")$cluster
    tr=Hybrid.Tree(data=mydata,cluster=cls, model="nbinom") ## bulild a tree structure for the resulting k clusters
    pdf(file=paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/clustering/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[m], "_clustering_nK_27.pdf", sep=""))
    plotHybrid.Tree(merge=tr, cluster=cls, logFC=mydata$logFC, tree.title=NULL, colorful = TRUE) ## plot the tree structure
    dev.off()
    condition.count.table.with.clustering = cbind(condition.count.table, cls)
    write.table(condition.count.table.with.clustering, file = paste("/home/ji72/RNAseq/mega_RNA-seq_Jul_2015_nextseq/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[m], "_nK_27_clustering_values.txt", sep=""), quote =F, row.names=T, col.names=T)
} #Setting nK as 27 (3 x 3 x 3)

#2. Draw the representation plots and heatmap ****on a local machine****
rm(list=ls(all=TRUE))
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
library(ggplot2); library(grid); library(gridExtra)
#Plot clustering results from model-based methods using the negative binomial model and K=27.
for (s in 1:length(list.of.bacteria.name)){
    bacteria.w.clustering.value = read.table(paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/Mega_RNA-seq_count_of_all_samples_Apr_2016_", list.of.bacteria.name[s], "_clustering_values.txt", sep=""), header=T)
    
    if (s == 2){ #for E.coli missing one of the 36hr replicates
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:5],1,mean),apply(bacteria.w.clustering.value[,6:8],1,mean),apply(bacteria.w.clustering.value[,9:11],1,mean))
        mean.bacteria.log10 = log10(mean.bacteria + 1)
        mean.bacteria.log10.w.value = as.data.frame(cbind(mean.bacteria.log10, bacteria.w.clustering.value[,12]))
    }
    else {
        mean.bacteria = cbind(apply(bacteria.w.clustering.value[,1:3],1,mean),apply(bacteria.w.clustering.value[,4:6],1,mean),apply(bacteria.w.clustering.value[,7:9],1,mean),apply(bacteria.w.clustering.value[,10:12],1,mean))
        mean.bacteria.log10 = log10(mean.bacteria + 1)
        mean.bacteria.log10.w.value = as.data.frame(cbind(mean.bacteria.log10, bacteria.w.clustering.value[,13]))
    }
    colnames(mean.bacteria.log10.w.value) = c("0hr","12hr","36hr","5.5d","cluster_number")
    zero = data.frame(mean.bacteria.log10.w.value[,1]); twelve = data.frame(mean.bacteria.log10.w.value[,2]); thirty.six = data.frame(mean.bacteria.log10.w.value[,3]); five.half = data.frame(mean.bacteria.log10.w.value[,4])
    zero.tag = cbind(zero, rep("0hr",dim(zero)[1]), mean.bacteria.log10.w.value[,5]); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]), mean.bacteria.log10.w.value[,5])
    ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1]), mean.bacteria.log10.w.value[,5]); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]), mean.bacteria.log10.w.value[,5])
    colnames(zero.tag) = c("log10Exp", "Time","Cluster");colnames(twelve.tag) = c("log10Exp", "Time","Cluster"); colnames(ts.tag) = c("log10Exp", "Time","Cluster"); colnames(fs.tag) = c("log10Exp", "Time","Cluster")
    total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.bacteria.log10.w.value),4)); total = cbind(total,gene.names)
    total$Time <- factor(total$Time); total$Cluster = factor(total$Cluster); total$gene.names = factor(total$gene.names)
    
    #Plot 1: First 15 clusters
    c1 = total[which(total$Cluster == 1),]; c2 = total[which(total$Cluster == 2),];c3 = total[which(total$Cluster == 3),]
    c4 = total[which(total$Cluster == 4),]; c5 = total[which(total$Cluster == 5),];c6 = total[which(total$Cluster == 6),]
    c7 = total[which(total$Cluster == 7),]; c8 = total[which(total$Cluster == 8),];c9 = total[which(total$Cluster == 9),]
    c10 = total[which(total$Cluster == 10),]; c11 = total[which(total$Cluster == 11),];c12 = total[which(total$Cluster == 12),]
    c13 = total[which(total$Cluster == 13),]; c14 = total[which(total$Cluster == 14),];c15 = total[which(total$Cluster == 15),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_first_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()
    
    #Plot 2: Biggest 15 clusters (This isn't very useful because there are too many entries)
    order.of.c.size = order(summary(total$Cluster),decreasing =T)
    c1 = total[which(total$Cluster == order.of.c.size[1]),]; c2 = total[which(total$Cluster == order.of.c.size[2]),];c3 = total[which(total$Cluster == order.of.c.size[3]),]
    c4 = total[which(total$Cluster == order.of.c.size[4]),]; c5 = total[which(total$Cluster == order.of.c.size[5]),];c6 = total[which(total$Cluster == order.of.c.size[6]),]
    c7 = total[which(total$Cluster == order.of.c.size[7]),]; c8 = total[which(total$Cluster == order.of.c.size[8]),];c9 = total[which(total$Cluster == order.of.c.size[9]),]
    c10 = total[which(total$Cluster == order.of.c.size[10]),]; c11 = total[which(total$Cluster == order.of.c.size[11]),]; c12 = total[which(total$Cluster == order.of.c.size[12]),]
    c13 = total[which(total$Cluster == order.of.c.size[13]),]; c14 = total[which(total$Cluster == order.of.c.size[14]),]; c15 = total[which(total$Cluster == order.of.c.size[15]),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()
    
    #     plot.list = list()
    #     for (p in 1:15){
    #         c = total[which(total$Cluster == order.of.c.size[p]),]
    #         plot.ind = ggplot(data =c, aes(x=c$Time, y = c$log10Exp, group = c$gene.names)) + geom_line() + theme_classic()
    #         plot.list[[p]] = plot.ind
    #     }
    #     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
    #     multiplot(plot.list[[1]], cols=5)
    #     dev.off()
    
}


#################################################################################################################
##############################################
#Option 5: EBSeq-HMM
##############################################
#################################################################################################################
#5/17/2016
#Done in my laptop
rm(list=ls(all=TRUE))

#Read in the mega table
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/")
CountTable <- read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt") #from deseq.for.samples.with.time.pca.for.mega.rna.R (a simple sorting >5 reads code that has nothing to do with deseq yielded this result)
CountTable = CountTable[,c(2:35, 37:70, 72:105)] 

#Label the samples accordingly
unchallenged.0hr = CountTable[,c(1:3,35:37,69:71)]; clean.prick.12hr = CountTable[,c(4,38,72)]; clean.prick.36hr = CountTable[,c(5,39,73)];clean.prick.5.5d = CountTable[,c(6,40,74)]; M.luteus.12hr = CountTable[,c(7,41,75)]; M.luteus.36hr = CountTable[,c(8,42,76)]; M.luteus.5.5d = CountTable[,c(9,43,77)]; E.coli.12hr = CountTable[,c(10,44,78)]; E.coli.36hr = CountTable[,c(11,79)]; 
E.coli.5.5d = CountTable[,c(12,46,80)]; S.mar.type.12hr = CountTable[,c(13,47,81)]; S.mar.type.36hr = CountTable[,c(14,48,82)]; S.mar.type.5.5d = CountTable[,c(15,49,83)]; E.fae.live.12hr = CountTable[,c(16,50,84)]; E.fae.live.36hr = CountTable[,c(17,51,85)]; E.fae.live.5.5d = CountTable[,c(18,52,86)]; P.rett.live.12hr = CountTable[,c(19,53,87)];
P.rett.live.36hr = CountTable[,c(20,54,88)]; P.rett.live.5.5d = CountTable[,c(21,55,89)]; Ecc15.12hr = CountTable[,c(22,56,90)]; Ecc15.36hr = CountTable[,c(23,57,91)]; Ecc15.5.5d = CountTable[,c(24,58,92)]; S.aureus.12hr = CountTable[,c(25,59,93)]; P.sneebia.12hr = CountTable[,c(26,60,94)]; S.mar.Db11.12hr = CountTable[,c(27,61,95)];
P.ento.12hr = CountTable[,c(28,62,96)]; E.fae.heatkilled.12hr = CountTable[,c(29,63,97)]; E.fae.heatkilled.36hr = CountTable[,c(30,64,98)]; E.fae.heatkilled.5.5d = CountTable[,c(31,65,99)]; P.rett.heatkilled.12hr = CountTable[,c(32,66,100)]; P.rett.heatkilled.36hr = CountTable[,c(33,67,101)]; P.rett.heatkilled.5.5d = CountTable[,c(34,68,102)]
unchallenged.average = cbind(apply(unchallenged.0hr[,1:3], 1, mean),apply(unchallenged.0hr[,4:6],1, mean),apply(unchallenged.0hr[,7:9], 1, mean))
colnames(unchallenged.average) = c("UC.Rep1", "UC.Rep2", "UC.Rep3")

M.luteus = cbind(unchallenged.average,M.luteus.12hr,M.luteus.36hr,M.luteus.5.5d) #M.luteus
E.coli = cbind(unchallenged.average,E.coli.12hr,E.coli.36hr,E.coli.5.5d) #E.coli
S.mar.type = cbind(unchallenged.average,S.mar.type.12hr,S.mar.type.36hr,S.mar.type.5.5d) #S.mar.type
E.fae.live = cbind(unchallenged.average,E.fae.live.12hr,E.fae.live.36hr,E.fae.live.5.5d) #E.fae.live
P.rett.live = cbind(unchallenged.average,P.rett.live.12hr,P.rett.live.36hr,P.rett.live.5.5d) #P.rett.live
Ecc15 = cbind(unchallenged.average,Ecc15.12hr,Ecc15.36hr,Ecc15.5.5d) #Ecc15

#list of conditions
list.of.conditions <- list(M.luteus, E.coli, S.mar.type, E.fae.live, P.rett.live,Ecc15)
list.of.name.of.conditions <- list("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
#de.of.all.conditions <- matrix(rownames(CountTable), nrow=dim(CountTable)[1], ncol=1)

#Execute EBSeq-HMM on 12hr, 36hr, and 5.5d time points. I excluded 0hr because it messes up the first point.
library("EBSeqHMM")
for (j in 1:length(list.of.conditions)){
    cat("EBSeq-HMM - We are working on the following comparison: 0hr-12hr-36hr-5.5d of ", as.character(list.of.name.of.conditions[j]),"samples \n")
    condition.count.table <- as.matrix(data.frame(list.of.conditions[j]))
    count.table <- condition.count.table
    Sizes <- MedianNorm(count.table)
    count.table <- GetNormalizedMat(count.table, Sizes)
    
    if (j == 2){ #getting rid of Rep2 of E.coli 36hr
        CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d")} #removing 36hr rep2 sample 
    else {CondVector <- c("t.0hr","t.0hr","t.0hr","t.12hr","t.12hr","t.12hr","t.36hr","t.36hr","t.36hr","t.5.5d","t.5.5d","t.5.5d") }    
    Conditions <- factor(CondVector, levels=c("t.0hr","t.12hr","t.36hr","t.5.5d"))
    
    EBSeqHMMGeneOut <- EBSeqHMMTest(Data=count.table, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=2)
    EBSeqHMMGeneOut$MgAllMAPChar #classifies the pattern of gene expression over time per gene
    AllPaths = GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=TRUE) #Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
    gene.list.for.a.given.path = GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, cutoff=0.8, OnlyDynamic=TRUE,Paths=NULL) #Obtain confident gene calls for classifying genes into expression paths, PP>=0.5 is most likely path
    number.of.genes.in.each.category = gene.list.for.a.given.path$NumEach; number.of.genes.in.each.category #already sorted
    genes.in.each.category = gene.list.for.a.given.path$EachPathNames; genes.in.each.category
    de.list = GetDECalls(EBSeqHMMGeneOut,FDR=.05) #Another way to obtain DE gene/isoform list at a certain FDR (but not PPs)
    
    #For a specific gene, draw a plot of expression
    #PlotExp(count.table, Conditions, "Gene_1")
    
    write.table(de.list, file=paste("EBSeqHMM_results_for_",list.of.name.of.conditions[j],".txt", sep=""), quote=F)
}

#Draw representation plot ****on a local machine****
list.of.bacteria.name = c("M.luteus", "E.coli", "S.mar.type", "E.fae.live", "P.rett.live", "Ecc15")
library(ggplot2); library(grid); library(gridExtra)
#Plot clustering results from EBSeq-HMM methods

for (s in 1:length(list.of.bacteria.name)){
    most.likely.paths = read.table(paste("clustering/EBSeq-HMM/EBSeqHMM_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)


    colnames(mean.bacteria.log10.w.value) = c("0hr","12hr","36hr","5.5d","cluster_number")
    zero = data.frame(mean.bacteria.log10.w.value[,1]); twelve = data.frame(mean.bacteria.log10.w.value[,2]); thirty.six = data.frame(mean.bacteria.log10.w.value[,3]); five.half = data.frame(mean.bacteria.log10.w.value[,4])
    zero.tag = cbind(zero, rep("0hr",dim(zero)[1]), mean.bacteria.log10.w.value[,5]); twelve.tag = cbind(twelve, rep("12hr",dim(twelve)[1]), mean.bacteria.log10.w.value[,5])
    ts.tag = cbind(thirty.six, rep("36hr",dim(thirty.six)[1]), mean.bacteria.log10.w.value[,5]); fs.tag = cbind(five.half, rep("5.5d",dim(five.half)[1]), mean.bacteria.log10.w.value[,5])
    colnames(zero.tag) = c("log10Exp", "Time","Cluster");colnames(twelve.tag) = c("log10Exp", "Time","Cluster"); colnames(ts.tag) = c("log10Exp", "Time","Cluster"); colnames(fs.tag) = c("log10Exp", "Time","Cluster")
    total = rbind(zero.tag, twelve.tag, ts.tag, fs.tag); gene.names=c(rep(row.names(mean.bacteria.log10.w.value),4)); total = cbind(total,gene.names)
    total$Time <- factor(total$Time); total$Cluster = factor(total$Cluster); total$gene.names = factor(total$gene.names)
    
    
    
    
    
    
    #Plot 1: First 15 clusters
    c1 = total[which(total$Cluster == 1),]; c2 = total[which(total$Cluster == 2),];c3 = total[which(total$Cluster == 3),]
    c4 = total[which(total$Cluster == 4),]; c5 = total[which(total$Cluster == 5),];c6 = total[which(total$Cluster == 6),]
    c7 = total[which(total$Cluster == 7),]; c8 = total[which(total$Cluster == 8),];c9 = total[which(total$Cluster == 9),]
    c10 = total[which(total$Cluster == 10),]; c11 = total[which(total$Cluster == 11),];c12 = total[which(total$Cluster == 12),]
    c13 = total[which(total$Cluster == 13),]; c14 = total[which(total$Cluster == 14),];c15 = total[which(total$Cluster == 15),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_first_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()
    
    #Plot 2: Biggest 15 clusters (This isn't very useful because there are too many entries)
    order.of.c.size = order(summary(total$Cluster),decreasing =T)
    c1 = total[which(total$Cluster == order.of.c.size[1]),]; c2 = total[which(total$Cluster == order.of.c.size[2]),];c3 = total[which(total$Cluster == order.of.c.size[3]),]
    c4 = total[which(total$Cluster == order.of.c.size[4]),]; c5 = total[which(total$Cluster == order.of.c.size[5]),];c6 = total[which(total$Cluster == order.of.c.size[6]),]
    c7 = total[which(total$Cluster == order.of.c.size[7]),]; c8 = total[which(total$Cluster == order.of.c.size[8]),];c9 = total[which(total$Cluster == order.of.c.size[9]),]
    c10 = total[which(total$Cluster == order.of.c.size[10]),]; c11 = total[which(total$Cluster == order.of.c.size[11]),]; c12 = total[which(total$Cluster == order.of.c.size[12]),]
    c13 = total[which(total$Cluster == order.of.c.size[13]),]; c14 = total[which(total$Cluster == order.of.c.size[14]),]; c15 = total[which(total$Cluster == order.of.c.size[15]),]
    plot1 = ggplot(data =c1, aes(x=c1$Time, y = c1$log10Exp, group = c1$gene.names)) + geom_line() + theme_classic()
    plot2 = ggplot(data =c2, aes(x=c2$Time, y = c2$log10Exp, group = c2$gene.names)) + geom_line() + theme_classic()
    plot3 = ggplot(data =c3, aes(x=c3$Time, y = c3$log10Exp, group = c3$gene.names)) + geom_line() + theme_classic()
    plot4 = ggplot(data =c4, aes(x=c4$Time, y = c4$log10Exp, group = c4$gene.names)) + geom_line() + theme_classic()
    plot5 = ggplot(data =c5, aes(x=c5$Time, y = c5$log10Exp, group = c5$gene.names)) + geom_line() + theme_classic()
    plot6 = ggplot(data =c6, aes(x=c6$Time, y = c6$log10Exp, group = c6$gene.names)) + geom_line() + theme_classic()
    plot7 = ggplot(data =c7, aes(x=c7$Time, y = c7$log10Exp, group = c7$gene.names)) + geom_line() + theme_classic()
    plot8 = ggplot(data =c8, aes(x=c8$Time, y = c8$log10Exp, group = c8$gene.names)) + geom_line() + theme_classic()
    plot9 = ggplot(data =c9, aes(x=c9$Time, y = c9$log10Exp, group = c9$gene.names)) + geom_line() + theme_classic()
    plot10 = ggplot(data =c10, aes(x=c10$Time, y = c10$log10Exp, group = c10$gene.names)) + geom_line() + theme_classic()
    plot11 = ggplot(data =c11, aes(x=c11$Time, y = c11$log10Exp, group = c11$gene.names)) + geom_line() + theme_classic()
    plot12 = ggplot(data =c12, aes(x=c12$Time, y = c12$log10Exp, group = c12$gene.names)) + geom_line() + theme_classic()
    plot13 = ggplot(data =c13, aes(x=c13$Time, y = c13$log10Exp, group = c13$gene.names)) + geom_line() + theme_classic()
    plot14 = ggplot(data =c14, aes(x=c14$Time, y = c14$log10Exp, group = c14$gene.names)) + geom_line() + theme_classic()
    plot15 = ggplot(data =c15, aes(x=c15$Time, y = c15$log10Exp, group = c15$gene.names)) + geom_line() + theme_classic()
    
    pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
    grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,plot15, ncol=5)
    dev.off()
    
    #     plot.list = list()
    #     for (p in 1:15){
    #         c = total[which(total$Cluster == order.of.c.size[p]),]
    #         plot.ind = ggplot(data =c, aes(x=c$Time, y = c$log10Exp, group = c$gene.names)) + geom_line() + theme_classic()
    #         plot.list[[p]] = plot.ind
    #     }
    #     pdf(file=paste("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.nk.50.model.based/",list.of.bacteria.name[s],"_MB_nK50_biggest_15_clustering_representation_plot.pdf", sep=""))
    #     multiplot(plot.list[[1]], cols=5)
    #     dev.off()
    
}

#For a given gene in the most_likely_path, identify the expression values in count.table and cbind
count.table = read.table("mega_RNA-seq_count_of_all_samples_Apr_2016.txt", header=T)
s=1
most.likely.paths = read.table(paste("EBSeqHMM_DE_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
#most.likely.paths = read.table(paste("clustering/EBSeq-HMM/EBSeqHMM_results_for_", list.of.bacteria.name[s], ".txt", sep=""), header=T)
exp.list = c()
for (n in 1:dim(most.likely.paths)[1]){
    gene.to.search = row.names(most.likely.paths[n,])
    gene.found = M.luteus[which(row.names(M.luteus) == gene.to.search),]
    exp.list = rbind(exp.list, gene.found)
}
mean.bacteria = cbind(apply(exp.list[,1:3],1,mean),apply(exp.list[,4:6],1,mean),apply(exp.list[,7:9],1,mean),apply(exp.list[,10:12],1,mean))
mean.bacteria.log10 = log10(mean.bacteria + 1)
most.likely.paths.M.luteus = cbind(most.likely.paths,mean.bacteria.log10)

par(mfrow=c(1,3))
most.likely.paths.M.luteus.c1 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "Up-Up-Up"),]
most.likely.paths.M.luteus.c1.exp.only = most.likely.paths.M.luteus.c1[,3:6]
colnames(most.likely.paths.M.luteus.c1.exp.only) = c("0hr","12hr","36hr","5.5d")
most.likely.paths.M.luteus.c1.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c1.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c1.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
ggplot(data = most.likely.paths.M.luteus.c1.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c1.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c1.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c1.exp.only.reshaped$id)) + geom_line() + theme_classic()

most.likely.paths.M.luteus.c2 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "Up-Up-Down"),]
most.likely.paths.M.luteus.c2.exp.only = most.likely.paths.M.luteus.c2[,3:6]
colnames(most.likely.paths.M.luteus.c2.exp.only) = c("0hr","12hr","36hr","5.5d")
most.likely.paths.M.luteus.c2.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c2.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c2.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
ggplot(data = most.likely.paths.M.luteus.c2.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c2.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c2.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c2.exp.only.reshaped$id)) + geom_line() + theme_classic()

most.likely.paths.M.luteus.c3 = most.likely.paths.M.luteus[which(most.likely.paths.M.luteus$Most_Likely_Path == "-Down-Down"),]
most.likely.paths.M.luteus.c3.exp.only = most.likely.paths.M.luteus.c3[,3:6]
colnames(most.likely.paths.M.luteus.c3.exp.only) = c("0hr","12hr","36hr","5.5d")
most.likely.paths.M.luteus.c3.exp.only.reshaped = reshape(most.likely.paths.M.luteus.c3.exp.only, direction="long", varying = 1:4, idvar = "id", ids = rownames(most.likely.paths.M.luteus.c3.exp.only), timevar="Time",v.names="log10.exp.value", times = c("0hr", "12hr", "36hr", "5.5d"))
ggplot(data = most.likely.paths.M.luteus.c3.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c3.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c3.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c3.exp.only.reshaped$id)) + geom_line() + theme_classic()


par(mfrow=c(1,3))
ggplot(data = most.likely.paths.M.luteus.c1.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c1.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c1.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c1.exp.only.reshaped$id)) + geom_line() + theme_classic()
ggplot(data = most.likely.paths.M.luteus.c2.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c2.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c2.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c2.exp.only.reshaped$id)) + geom_line() + theme_classic()
ggplot(data = most.likely.paths.M.luteus.c3.exp.only.reshaped , aes(x=most.likely.paths.M.luteus.c3.exp.only.reshaped$Time, y = most.likely.paths.M.luteus.c3.exp.only.reshaped$log10.exp.value, group = most.likely.paths.M.luteus.c3.exp.only.reshaped$id)) + geom_line() + theme_classic()


#################################################################################################################
##############################################
#Option 6: Do the clustering on 12hr, 36hr, and 5.5d samples excluding the 0hr samples
##############################################
#################################################################################################################

#In this case, just 9 clusters