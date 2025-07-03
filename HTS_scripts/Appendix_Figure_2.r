#
# R syntax to reproduce information for Appendix Figure 2 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al.
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

sero <- read.csv(paste(dir, "paired_sera.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

# pH1N1
pcr.ph1 <- unique(hc[hc$swine==1&!is.na(hc$swine),1:2])
sero.ph1 <- merge(sero[c(1:4,8,9)],pcr.ph1,by=c("hhID","member"))
sero.ph1$neut.ph1.paired <- sero.ph1$neut_pH1_c/sero.ph1$neut_pH1
sero.ph1$hi.ph1.paired <- sero.ph1$HI_pH1_c/sero.ph1$HI_pH1
plotdata.ph1 <- data.frame(rise=c(1,2,4,8,16,32,64,128,256))
for(i in 1:9){
   plotdata.ph1$prop.neut[i] <- sum(log2(sero.ph1$neut.ph1.paired)+1>=i)/dim(sero.ph1)[1]
   plotdata.ph1$prop.hi[i] <- sum(log2(sero.ph1$hi.ph1.paired)+1>=i)/dim(sero.ph1)[1]
}

# sH1N1
pcr.sh1 <- unique(hc[hc$AH1==1&!is.na(hc$AH1),1:2])
sero.sh1 <- merge(sero[c(1,2,5,10)],pcr.sh1,by=c("hhID","member"))
sero.sh1$hi.sh1.paired <- sero.sh1$HI_sH1_c/sero.sh1$HI_sH1
plotdata.sh1 <- data.frame(rise=c(1,2,4,8,16,32,64,128,256))
for(i in 1:9){
   plotdata.sh1$prop.hi[i] <- sum(log2(sero.sh1$hi.sh1.paired)+1>=i)/dim(sero.sh1)[1]
}

# sH3N2
pcr.sh3 <- unique(hc[hc$AH3==1&!is.na(hc$AH3),1:2])
sero.sh3 <- merge(sero[c(1,2,6,7,11,12)],pcr.sh3,by=c("hhID","member"))
sero.sh3$neut.sh3.paired <- sero.sh3$neut_sH3_c/sero.sh3$neut_sH3
sero.sh3$hi.sh3.paired <- sero.sh3$HI_sH3_c/sero.sh3$HI_sH3
plotdata.sh3 <- data.frame(rise=c(1,2,4,8,16,32,64,128,256))
for(i in 1:9){
   plotdata.sh3$prop.neut[i] <- sum(log2(sero.sh3$neut.sh3.paired)+1>=i)/dim(sero.sh3)[1]
   plotdata.sh3$prop.hi[i] <- sum(log2(sero.sh3$hi.sh3.paired)+1>=i)/dim(sero.sh3)[1]
}

# 3-panel plot (pH1N1,sH1N1,sH3N2)

windows(width=6,height=12)
layout(matrix(1:3,ncol=1))

par(mar=c(4.5,4,1.1,1),xaxs="i",yaxs="i")
plot(NA,xlim=c(-0.2,8.1),ylim=c(0,1.05),axes=FALSE,xlab="", main="Pandemic H1N1",
     ylab="Proportion with equal or greater rise")
lines(log2(plotdata.ph1$rise),plotdata.ph1$prop.neut,lty=2)
lines(log2(plotdata.ph1$rise),plotdata.ph1$prop.hi)
axis(1,at=0:8,labels=c(expression(""<="1"),2,4,8,16,32,64,128,256))
axis(2,pos=0,at=0:5/5,labels=c("0%","20%","40%","60%","80%","100%"),las=1)
legend(5.4,0.8,legend=c("VN  ","HAI"),lty=2:1,cex=0.9)

par(mar=c(4.5,4,1.1,1),xaxs="i",yaxs="i")
plot(NA,xlim=c(-0.2,8.1),ylim=c(0,1.05),axes=FALSE,xlab="", main="Seasonal H1N1",
     ylab="Proportion with equal or greater rise")
lines(log2(plotdata.sh1$rise),plotdata.sh1$prop.hi)
axis(1,at=0:8,labels=c(expression(""<="1"),2,4,8,16,32,64,128,256))
axis(2,pos=0,at=0:5/5,labels=c("0%","20%","40%","60%","80%","100%"),las=1)

par(mar=c(4.5,4,1.1,1),xaxs="i",yaxs="i")
plot(NA,xlim=c(-0.2,8.1),ylim=c(0,1.05),axes=FALSE,xlab="Rise in antibody titer from baseline to convalescence", main="Seasonal H3N2",
     ylab="Proportion with equal or greater rise")
lines(log2(plotdata.sh3$rise),plotdata.sh3$prop.neut,lty=2)
lines(log2(plotdata.sh3$rise),plotdata.sh3$prop.hi,lty=1)
axis(1,at=0:8,labels=c(expression(""<="1"),2,4,8,16,32,64,128,256))
axis(2,pos=0,at=0:5/5,labels=c("0%","20%","40%","60%","80%","100%"),las=1)
legend(5.4,0.8,legend=c("VN  ","HAI"),lty=2:1,cex=0.9)

# End of script
