#Drawing a curve for fly survival (Acp65Aa RNAi)
#July 28th, 2016
#Joo Hyun Im
#exp: date of when the experiment started (trial #)
#day: the day that flies died
#status: 1 is dead, 0 is not dead
################################################################################

rm(list=ls(all=TRUE))
library(survival)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/RNAi_validation/")
data.snb<-read.table("Acp65Aa_RNAi_survival_against_P.sneebia_0.01_29C.txt", header=TRUE, sep="\t", as.is=c(factor,factor,numeric,numeric,factor))
data.rett<-read.table("Acp65Aa_RNAi_survival_against_P.rettgeri_0.1_29C.txt", header=TRUE, sep="\t", as.is=c(factor,factor,numeric,numeric,factor))

#Graph 1: All replicates combined
#P.sneebia
data.snb.gd = data.snb[grep("GD", data.snb$tag),]; data.snb.kk = data.snb[grep("KK", data.snb$tag),]
par(mfrow = c(1,2))
plot(survfit(Surv(data.snb.gd$day,data.snb.gd$status)~data.snb.gd$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.sneebia infection (OD 0.01)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("Acp65Aa RNAi (n=82)","control-GD (n=87)"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)
plot(survfit(Surv(data.snb.kk$day,data.snb.kk$status)~data.snb.kk$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.sneebia infection (OD 0.01)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("Acp65Aa RNAi (n=81)","control-KK (n=74)"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)

model.snb.gd<-coxph(Surv(data.snb.gd$day,data.snb.gd$status)~data.snb.gd$tag); anova(model.snb.gd)
model.snb.kk<-coxph(Surv(data.snb.kk$day,data.snb.kk$status)~data.snb.kk$tag); anova(model.snb.kk)

#P.rett
data.rett.gd = data.rett[grep("GD", data.rett$tag),]; data.rett.kk = data.rett[grep("KK", data.rett$tag),]
par(mfrow = c(1,2))
plot(survfit(Surv(data.rett.gd$day,data.rett.gd$status)~data.rett.gd$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.rettgeri infection (OD 0.1)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("Acp65Aa RNAi (n=41)","control-GD (n=43)"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)
plot(survfit(Surv(data.rett.kk$day,data.rett.kk$status)~data.rett.kk$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.rettgeri infection (OD 0.1)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("Acp65Aa RNAi (n=43)","control-KK (n=40)"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)

model.rett.gd<-coxph(Surv(data.rett.gd$day,data.rett.gd$status)~data.rett.gd$tag); anova(model.rett.gd)
model.rett.kk<-coxph(Surv(data.rett.kk$day,data.rett.kk$status)~data.rett.kk$tag); anova(model.rett.kk)


#################################

#Graph 3: All replicates separated
par(mfrow = c(1,2))
plot(survfit(Surv(data.snb$day,data.snb$status)~data.snb$tag+data.snb$rep),col=c("darkred", "darkred","darkred","darkred","darkred","darkred","grey50","grey50","grey50","grey50","grey50","grey50"),lty=c(1,2,3,4,5,6,1,2,3,4,5,6),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", main="P.sneebia infection (OD 0.01)", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("1","2","3","4","5","6"),lty=c(1,2,3,4,5,6), lwd=3, cex=1.2)
plot(survfit(Surv(data.rett$day,data.rett$status)~data.rett$tag+data.rett$rep),col=c("darkred", "darkred","darkred","darkred","grey50","grey50","grey50","grey50"),lty=c(1,2,3,4,1,2,3,4),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", main="P.rettgeri infection (OD 1)", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("1","2","3","4"),lty=c(1,2,3,4), lwd=3, cex=1.2)


#Graph 4: All replicates combined -- exclude one replicate of P.rettgeri that had a lift in survival (rep 1)
data.rett2 = data.rett[which(data.rett$rep !=1),]
par(mfrow = c(1,2))
plot(survfit(Surv(data.snb$day,data.snb$status)~data.snb$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.sneebia infection (OD 0.01)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("NucB1 RNAi","Driver control"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)
plot(survfit(Surv(data.rett2$day,data.rett2$status)~data.rett2$tag),col=c("darkred", "grey50"),lty=c(1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.rettgeri infection (OD 1)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("NucB1 RNAi","Driver control"),col=c("darkred", "grey50"), lty=c(1,1), lwd=3, cex=1.2)

model.snb<-coxph(Surv(data.snb$day,data.snb$status)~data.snb$tag); anova(model.snb)
model.rett2<-coxph(Surv(data.rett2$day,data.rett2$status)~data.rett2$tag); anova(model.rett2)


#Graph 3: by condition (all reps separated)
par(mfrow = c(1,3))
data.rnai.only = data[which(data$tag == "Nucb1_RNAi"),]
plot(survfit(Surv(data.rnai.only$day,data.rnai.only$status)~data.rnai.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("NucB1 RNAi 1","NucB1 RNAi 2","NucB1 RNAi 3","NucB1 RNAi 4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)
data.cont1.only = data[which(data$tag == "Nucb1_RNAi_bg1"),]
plot(survfit(Surv(data.cont1.only$day,data.cont1.only$status)~data.cont1.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("Control1-1","Control1-2","Control1-3","Control1-4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)
data.cont2.only = data[which(data$tag == "Nucb1_RNAi_bg2"),]
plot(survfit(Surv(data.cont2.only$day,data.cont2.only$status)~data.cont2.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("Control2-1","Control2-2","Control2-3","Control2-4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)


#Graph 6: by condition (all reps separated)
par(mfrow = c(1,3))
data.rnai.only = data[which(data$tag == "Nucb1_RNAi"),]
plot(survfit(Surv(data.rnai.only$day,data.rnai.only$status)~data.rnai.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("NucB1 RNAi 1","NucB1 RNAi 2","NucB1 RNAi 3","NucB1 RNAi 4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)
data.cont1.only = data[which(data$tag == "Nucb1_RNAi_bg1"),]
plot(survfit(Surv(data.cont1.only$day,data.cont1.only$status)~data.cont1.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("Control1-1","Control1-2","Control1-3","Control1-4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)
data.cont2.only = data[which(data$tag == "Nucb1_RNAi_bg2"),]
plot(survfit(Surv(data.cont2.only$day,data.cont2.only$status)~data.cont2.only$rep),col=c("darkred", "blue","chartreuse","burlywood4"),lty=c(1,1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5)
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft",inset=0.05, bty="n", c("Control2-1","Control2-2","Control2-3","Control2-4"),col=c("darkred", "blue","chartreuse","burlywood4"), lty=c(1,1,1,1), lwd=3, cex=1.2)


#Graph 7: P.sneebia and P.rettgeri infection together: All replicates combined
rm(list=ls(all=TRUE))
library(survival)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/RNAi_validation/")
data.snb<-read.table("RNAi_survival_against_P.sneebia_0.01_29C.txt", header=TRUE, sep="\t",as.is=c(factor,factor,numeric,numeric,factor))
data.rett<-read.table("RNAi_survival_against_P.rettgeri_1_29C.txt", header=TRUE, sep="\t",as.is=c(factor,factor,numeric,numeric,factor))

par(mfrow = c(1,2))
plot(survfit(Surv(data.snb$day,data.snb$status)~data.snb$tag),col=c("darkred", "grey50", "grey70"),lty=c(1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.sneebia infection (OD 0.01)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("NucB1 RNAi (n=86)","Driver control (n=76)","RNAi control (n=79)"),col=c("darkred", "grey50", "grey70"), lty=c(1,1,1), lwd=3, cex=1.2)
plot(survfit(Surv(data.rett$day,data.rett$status)~data.rett$tag),col=c("darkred", "grey50", "grey70"),lty=c(1,1,1),lwd=4,cex=1.7,yaxt='n',xaxt='n',cex.main=1.5,cex.axis=1.5, ylab="Proportion Alive", xlab="Days Post Infection", cex.lab=1.5, main ="P.rettgeri infection (OD 1)")
axis(2,las=2,cex.axis=1.4,lwd=3)
axis(1,lwd=3, cex.axis=1.4)
box(lwd=3)
legend("bottomleft", bty="n", c("NucB1 RNAi (n=82)","Driver control (n=81)","RNAi control (n=84)"),col=c("darkred", "grey50", "grey70"), lty=c(1,1,1), lwd=3, cex=1.2)

#Graph 8: P.sneebia and P.rettgeri infection together: All replicates combined, RNAi control (bg2) is dropped
rm(list=ls(all=TRUE))
library(survival)
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/RNAi_validation/")
data.snb<-read.table("RNAi_survival_against_P.sneebia_0.01_29C.txt", header=TRUE, sep="\t",as.is=c(factor,factor,numeric,numeric,factor))
data.rett<-read.table("RNAi_survival_against_P.rettgeri_1_29C.txt", header=TRUE, sep="\t",as.is=c(factor,factor,numeric,numeric,factor))
data.snb = data.snb[which(data.snb$tag != "Nucb1_RNAi_bg2"),]; data.rett = data.rett[which(data.rett$tag != "Nucb1_RNAi_bg2"),]
data.snb$tag = factor(data.snb$tag); data.rett$tag = factor(data.rett$tag)


