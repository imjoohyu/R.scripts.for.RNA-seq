#Drawing a curve for fly survival (Attp40 RNAi)
#May 19th, 2016
#Joo Hyun Im
#exp: sample#
#day: the day that flies died
#status: 1 is dead, 0 is not dead
################################################################################
#All (5/19/2016)
rm(list=ls(all=TRUE))
library(survival)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/RNAi_validation/")
data<-read.table("Attp40_RNAi_survival_against_P.sneebia_all_data.txt", header=TRUE, sep="\t",as.is=c(factor,factor,factor,numeric,numeric,factor))

#Question 1: Which concentration of P. sneebia is suitable for infection?
#Plot PBS, 29C OD 0.01, 29C OD 0.1, and 29C OD 1
data.q1 = data[c(1:165),]; attach(data.q1); rep<-as.factor(rep); concentration<-as.factor(concentration)
plot(survfit(Surv(day,status)~concentration),col=c("grey","pink","brown1","darkred"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("PBS","29C OD 0.01","29C OD 0.1","29C OD 1"),col=c("grey","pink","brown1","darkred"), lty=c(1,1,1,1), lwd=3, cex=1.2) #a plot for fun

#Added more concentrations. Added on 6/9/2016:
data.multi.only <-read.table("Attp40_RNAi_survival_against_P.sneebia_multiple_conc_29C.txt", header=T, sep="\t",as.is=c(factor,factor,factor,numeric,numeric,factor))
attach(data.multi.only); concentration<-as.factor(concentration)
plot(survfit(Surv(day,status)~concentration),col=c("lightpink3","indianred2","firebrick3","indianred4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("29C OD 0.1","29C OD 0.5","29C OD 0.6", "29C OD 0.75"),col=c("lightpink3","indianred2","firebrick3","indianred4"), lty=c(1,1,1,1), lwd=3, cex=1.2) #a plot for fun


#Question 2: At OD 1, is there a difference between flies at 29C and flies at 25C?
#Plot consolidated 29C OD 1 and 25C OD 1
data.q2 = data[c(44:224),]; attach(data.q2); data.q2$temp<-as.factor(data.q2$temp)
plot(survfit(Surv(day,status)~temp),col=c("darkred","blue"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("topright",inset=0.05, bty="n" , c("29C OD 1","25C OD 1"),col=c("darkred","blue"),lty=c(1,1),lwd=3, cex=1.5) #a plot for fun

#Question 3: How noisy are the graphs at 29C?
#Plot separate 29C OD 1 graphs (age-controlled, date differs)
data.q3 = data[c(44:165),]; attach(data.q3); data.q3$exp<-as.factor(data.q3$exp)
plot(survfit(Surv(day,status)~exp),col=c("lightpink","lightcoral","lightcoral","indianred4","indianred4","indianred4","indianred4"),lty=c(1,1,1,1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3); axis(1,lwd=3, cex.axis=1.4); box(lwd=3)
legend("topright",inset=0.05, bty="n" , c("D1","D2","D2","D3","D3","D3","D3"),col=c("lightpink","lightcoral","lightcoral","indianred4","indianred4","indianred4","indianred4"),lty=c(1,1,1,1,1,1,1),lwd=3, cex=1.5) #a plot for fun

#Question 4: How noisy are the graphs at 25C?
#Plot separate 25C OD 1 graphs (date-controlled, age differs)
data.q4 = data[c(166:224),]; attach(data.q4); data.q4$exp<-as.factor(data.q4$exp)
plot(survfit(Surv(day,status)~exp),col=c("blue","blue","blue"),lty=c(1,2,3),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3); axis(1,lwd=3, cex.axis=1.4); box(lwd=3)
legend("topright",inset=0.05, bty="n", c("D1","D1","D1"), col=c("blue","blue","blue"),lty=c(1,2,3),lwd=3, cex=1.5) #a plot for fun

bacteria.load.table = matrix(log10(c(60,1280,540,20,20,1)))
#bacteria.load.table = cbind(bacteria.load.table, c("sample1", "sample2", "sample3","sample1", "sample2", "sample3"))
bacteria.load.table = cbind(bacteria.load.table, c("25C","25C","25C","29C","29C","29C"))
colnames(bacteria.load.table) = c("concentration","temperature")
bacteria.load.table = data.frame(bacteria.load.table); attach(bacteria.load.table)
library(ggplot2)
ggplot(data = bacteria.load.table, aes(x = bacteria.load.table$temperature, y = bacteria.load.table$concentration) + geom_bar(color=c("darkred"), width=0.8, stat="identity") + xlab("Sample") + ylab("Count/mL") + ggtitle("Amount of bacteria present in the survivors")
