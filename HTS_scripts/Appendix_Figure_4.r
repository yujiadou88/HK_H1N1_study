#
# R syntax to reproduce information for Appendix Figure 4 from:
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


swine <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
sh1 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH1==1])
sh3 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH3==1])


pandemic <- sero[sero$hhID%in%swine,c(1:4,8,9)]

# index with H1N1 pandemic
pand_index <- pandemic[pandemic$member==0,]
pand_index <- merge(pand_index,demog[c(1,2,5)],by=c("hhID","member"),all.x=TRUE)
pand_index$rise.neut <- pand_index$neut_pH1_c/pand_index$neut_pH1
pand_index$rise.hi <- pand_index$HI_pH1_c/pand_index$HI_pH1

# index with sH3N2
sH3N2_index <- sero[sero$hhID%in%sh3&sero$member==0,c(1,2,6,7,11,12)]
sH3N2_index <- merge(sH3N2_index,demog[c(1,2,5)],by=c("hhID","member"),all.x=TRUE)
sH3N2_index$rise.neut <- sH3N2_index$neut_sH3_c/sH3N2_index$neut_sH3
sH3N2_index$rise.hi <- sH3N2_index$HI_sH3_c/sH3N2_index$HI_sH3

#
# plot
#

windows(width=10,height=10)
layout(matrix(1:4,ncol=2,byrow=F))

set.seed(12345)
par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Pandemic H1N1")
points(jitter(pand_index$tamiflu+0.5,factor=0.5),jitter(log2(pand_index$rise.hi),factor=0.5),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)
mtext("GMT increase",side=2,line=3)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Pandemic H1N1")
points(jitter(pand_index$tamiflu+0.5,factor=0.5),jitter(log2(pand_index$rise.neut),factor=0.5),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)
mtext("GMT increase",side=2,line=3)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Seasonal H3N2")
points(jitter(sH3N2_index$tamiflu+0.5,factor=0.5),jitter(log2(sH3N2_index$rise.hi),factor=0.5),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Seasonal H3N2")
points(jitter(sH3N2_index$tamiflu+0.5,factor=0.5),jitter(log2(sH3N2_index$rise.neut),factor=0.5),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)

#
# order the points manually
#

windows(width=10,height=10)
layout(matrix(1:4,ncol=2,byrow=F))

set.seed(12345)
par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Pandemic H1N1")
points(rep(0.5,5),log2(c(8,16,32,64,128)),pch=16)
points(rep(1.5,3),log2(c(2,8,64)),pch=16)
points(c(rep(1.47,2),rep(1.53,2)),log2(c(1,4,1,4)),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)
mtext("GMT increase",side=2,line=3)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Pandemic H1N1")
points(rep(0.5,2),log2(c(32,64)),pch=16)
points(rep(1.5,3),log2(c(1,8,32)),pch=16)
points(c(0.45,0.5,0.55),log2(c(16,16,16)),pch=16)
points(c(rep(1.47,2),rep(1.53,2)),log2(c(4,16,4,16)),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)
mtext("GMT increase",side=2,line=3)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Seasonal H3N2")
points(rep(0.5,3),log2(c(2,4,16)),pch=16)
points(rep(1.5,6),log2(c(2,4,8,16,128,256)),pch=16)
points(c(0.47,0.53),log2(c(8,8)),pch=16)
points(c(0.45,0.5,0.55),log2(rep(1,3)),pch=16)
points(c(1.4,1.45,1.5,1.55,1.6),log2(rep(1,5)),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)

par(mar=c(4,4,1,1))
plot(NA,xlim=c(0,2),ylim=c(0,8),axes=FALSE,xlab="",ylab="",main="Seasonal H3N2")
points(0.5,log2(4),pch=16)
points(c(rep(0.47,2),rep(0.53,2)),log2(c(2,16,2,16)),pch=16)
points(c(rep(0.45,1),rep(0.5,1),rep(0.55,1)),log2(c(8,8,8)),pch=16)
points(rep(1.5,2),log2(c(1,8)),pch=16)
points(c(rep(1.47,3),rep(1.53,3)),log2(c(2,16,64,2,16,64)),pch=16)
points(c(rep(1.45,1),rep(1.5,1),rep(1.55,1)),log2(c(4,4,4)),pch=16)
axis(1,at=0:3-0.5,labels=c("","No oseltamivir","Oseltamivir",""))
axis(2,at=log2(c(1,2,4,8,16,32,64,128,256)),labels=c(1,2,4,8,16,32,64,128,256),las=1)

# End of script
