#Making figures for lab meeting 11/3/2015
#October 29th, 2015
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry

#1. Count the number of DEGs per infection-time condition
#read the file that has a number of DEGs per infection-time condition
num.degs = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/num.of.DEGs.Nov.2015.txt", header=F)
#Number from list_of_DEGs folder created by GO.analysis.pulling.out.sig.genes.
# tag = c(rep("DEG-up",31),rep("DEG-down",31)); 
condition = c(rep(c("clean.prick.12hr","M.luteus.12hr","E.coli.12hr","S.mar.type.12hr","Ecc15.12hr","P.rett.live.12hr","P.rett.heatkilled.12hr",
                                                             "E.fae.live.12hr", "E.fae.heatkilled.12hr","S.aureus.12hr", "P.sneebia.12hr", "S.mar.Db11.12hr", "P.ento.12hr","clean.prick.36hr",
                                                             "M.luteus.36hr","E.coli.36hr","S.mar.type.36hr","Ecc15.36hr","P.rett.live.36hr","P.rett.heatkilled.36hr",
                                                             "E.fae.live.36hr", "E.fae.heatkilled.36hr", "clean.prick.5.5d", "M.luteus.5.5d",   "E.coli.5.5d",
     "S.mar.type.5.5d",  "Ecc15.5.5d","P.rett.live.5.5d","P.rett.heatkilled.5.5d", "E.fae.live.5.5d", "E.fae.heatkilled.5.5d"),2))
# num.degs.redone = as.data.frame(cbind(condition, tag, num.degs[,2])); colnames(num.degs.redone) = c("condition","tag","number.of.genes")
# #data <- tapply(num.degs.redone$number.of.genes, list(num.degs.redone$tag,num.degs.redone$condition))
# #barplot(data,beside=T,main="Number of DEGs per infection conditions",xlab="Infection conditions",ylab="Number of genes")
num.degs.down = num.degs[1:31,]; num.degs.up = num.degs[32:62,]
par(mfrow=c(2,3)) 
barplot(num.degs.up[1:13,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), ylim=c(0,700), xlab="12hr",cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.up[14:22,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.up[23:31,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot
barplot(num.degs.down[1:13,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), ylim=c(0,700), xlab="12hr", cex.lab=1.5,cex.axis=1.5) #12hr barplot
barplot(num.degs.down[14:22,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="36hr",cex.lab=1.5,cex.axis=1.5) #36hr barplot
barplot(num.degs.down[23:31,2], col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1"),ylim=c(0,700), xlab="5.5day",cex.lab=1.5,cex.axis=1.5) #5.5d barplot
#Just get the legend
barplot(c(1:13), col=c("darkgoldenrod1","pink1","lightblue1","lightskyblue","royalblue1","deepskyblue3","deepskyblue2","firebrick2","coral1","coral4","dodgerblue3","dodgerblue4","darkslateblue"), 
        legend = c("clean.prick","M.luteus","E.coli","S.mar.type","Ecc15","P.rett.live","P.rett.heatkilled",
                   "E.fae.live", "E.fae.heatkilled","S.aureus", "P.sneebia", "S.mar.Db11", "P.ento"))
# ecc15 "royalblue1",ecoli "lightblue1",efae.live "firebrick2",efae.hk "coral1",mlu "pink1",p.ento "darkslateblue",
# p.rett.live "deepskyblue3",p.rett.hk "deepskyblue2", p.sneeb "dodgerblue3",s.aureus "coral4",smar.db11 "dodgerblue4",smartype"lightskyblue")


#2. Counting the frequency of GO terms in DAVID_GO_GEA folder using UNIX
#cd Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_upregulated_DE/
#cat GO_term_only_*.txt > GO_term_only_upregulated_all_conditions.txt
#cut -f1 -d'   ' GO_term_only_upregulated_all_conditions.txt > GO_term_only_upregulated_all_conditions_GO_only.txt
#sort GO_term_only_upregulated_all_conditions_GO_only.txt | uniq -c | sort -nr > GO_term_only_upregulated_all_conditions_GO_only_uniq.txt #control + v + and then tab
#Get rid of term 30
#sort GO_term_only_upregulated_all_conditions_GO_only.txt | uniq | cut -f1 -d "~" > GO_term_only_upregulated_all_conditions_GO_only_uniq_no_count.txt
#cut -f1 -d ":" GO_term_only_upregulated_all_conditions_GO_only_uniq.txt > GO_term_only_upregulated_all_conditions_GO_only_uniq_count_only.txt
#Open the GO_term_only_upregulated_all_conditions_GO_only_uniq_count_only.txt file in Excel and delete the second column.

#change up/down first!
#go.terms.counts = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_upregulated_DE/GO_term_only_upregulated_all_conditions_GO_only_uniq_count_only.txt", header=F)
go.terms.counts = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_downregulated_DE/GO_term_only_downregulated_all_conditions_GO_only_uniq_count_only.txt", header=F)
library('plyr'); count(go.terms.counts)
int.hist = function(x,ylab="Frequency",...) {
    barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1) }
int.hist(as.matrix(go.terms.counts), xlab="Number of infection conditions", 
         ylab="Number of GO terms", xlim=c(0,31), ylim=c(0,80), cex.lab=1.5, cex.axis=1.5, cex.main=1.5) #main = "Frequency of up DEG GO terms across infection conditions", 

#3. Counting the frequency of KEGG terms in DAVID_GO_GEA folder using UNIX
#cd Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/DAVID_GO_GEA/GO_upregulated_DE/
#cat KEGG_GEA_term_only_*.txt > KEGG_GEA_term_only_upregulated_all_conditions.txt

#cut -f1 -d'   ' KEGG_GEA_term_only_upregulated_all_conditions.txt > KEGG_GEA_term_only_upregulated_all_conditions_GO_only.txt
#sort KEGG_GEA_term_only_upregulated_all_conditions_GO_only.txt | uniq -c | sort -nr > KEGG_GEA_term_only_upregulated_all_conditions_GO_only_uniq.txt
#sort KEGG_GEA_term_only_upregulated_all_conditions_GO_only.txt | uniq | cut -f1 -d "~" > KEGG_GEA_term_only_upregulated_all_conditions_GO_only_uniq_no_count.txt



