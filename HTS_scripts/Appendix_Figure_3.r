#
# R syntax to reproduce information for Appendix Figure 3 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al.
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

sero <- read.csv(paste(dir, "paired_sera.csv", sep=""))
demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

# Define infection (PCR+)
hc2 <- reshape(hc[1:4],timevar="visit",idvar=c("hhID","member"),direction="wide")
hc.contact <- hc2[hc2$member>0,]
hc.contact$infect <- 0
hc.contact$infect[hc.contact$qPCR.1==1|(!is.na(hc.contact$qPCR.2)&hc.contact$qPCR.2>0)
                      |(!is.na(hc.contact$qPCR.3)&hc.contact$qPCR.3>0)] <- 1
serohc <- merge(hc.contact[c(1,2,6)],sero[c(1:3,5,6)],by=c("hhID","member"))
serohc <- merge(serohc,demog[1:3],by=c("hhID","member"),all.x=TRUE)

ph1 <- unique(hc$hhID[hc$swine==1])
sh1 <- unique(hc$hhID[hc$AH1==1])
sh3 <- unique(hc$hhID[hc$AH3==1])

plot.ph1 <- serohc[serohc$hhID%in%ph1,]
plot.sh1 <- serohc[serohc$hhID%in%sh1,]
plot.sh3 <- serohc[serohc$hhID%in%sh3,]
plot.ph1 <- plot.ph1[order(plot.ph1$infect),]
plot.sh1 <- plot.sh1[order(plot.sh1$infect),]
plot.sh3 <- plot.sh3[order(plot.sh3$infect),]

pch1 <- c(1,17)

windows(width=6,height=10)
layout(matrix(1:3,ncol=1))
par(mar=c(4.5,4,1,1))

set.seed(123456)
plot(NA,xlim=c(0,80),ylim=c(2,10.5),axes=FALSE,xlab="",ylab="Baseline antibody titer",main="")
points(jitter(plot.ph1$age),jitter(log2(plot.ph1$neut_pH1)),pch=pch1[plot.ph1$infect+1],col=plot.ph1$infect+1,cex=1.2)
axis(1,at=0:8*10)
axis(2,at=log2(c(5,10,20,40,80,160,320,640,1280)),labels=c("<1:10",10,20,40,80,160,320,640,1280),las=1)
legend(60,9.5,legend=c("Not infected","Infected"),col=1:2,pch=pch1)
mtext("Age",side=1,line=2.5,cex=0.7)
mtext("Pandemic H1N1",side=3,line=-1,font=2)

plot(NA,xlim=c(0,80),ylim=c(2,10.5),axes=FALSE,xlab="",ylab="Baseline antibody titer",main="")
points(jitter(plot.sh1$age),jitter(log2(plot.sh1$HI_sH1)),pch=pch1[plot.sh1$infect+1],col=plot.sh1$infect+1,cex=1.2)
axis(1,at=0:8*10)
axis(2,at=log2(c(5,10,20,40,80,160,320,640,1280)),labels=c("<1:10",10,20,40,80,160,320,640,1280),las=1)
mtext("Age",side=1,line=2.5,cex=0.7)
mtext("Seasonal H1N1",side=3,line=-1,font=2)

plot(NA,xlim=c(0,80),ylim=c(2,10.5),axes=FALSE,xlab="",ylab="Baseline antibody titer",main="")
points(jitter(plot.sh3$age),jitter(log2(plot.sh3$neut_sH3)),pch=pch1[plot.sh3$infect+1],col=plot.sh3$infect+1,cex=1.2)
axis(1,at=0:8*10)
axis(2,at=log2(c(5,10,20,40,80,160,320,640,1280)),labels=c("<1:10",10,20,40,80,160,320,640,1280),las=1)
mtext("Age",side=1,line=2.5,cex=0.7)
mtext("Seasonal H3N2",side=3,line=-1,font=2)

# End of script
