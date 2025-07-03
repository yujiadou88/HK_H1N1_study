#
# Script to reproduce information for Appendix Figure 1 from:
#
# Cowling BJ, Fang VJ, Riley S, Peiris JS, Leung GM. 
# An estimate of the serial interval of influenza using 
# laboratory-confirmed natural infections in households. 
# Epidemiology, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Feburary 2, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)

###---------------------------------------------------- Lab-confirmed secondary cases ------------------------------------------------###

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
for (i in 1:nrow(hc)){
     if(!is.na(hc$PCR[i]) & (is.na(hc$culture[i])|(!is.na(hc$culture[i])&hc$culture[i]==0)) ) hc$culture[i] <- hc$PCR[i]
}
hc <- hc[-5]

hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="culture")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V0","V1","V2","V3","V4")

for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
                  | (is.na(hculture$V1[i]) & hculture$V0[i]==0) | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
              {hculture$exd_index[i]<-1}    else {hculture$exd_index[i]<-0}
     if(hculture$member[i]==0 &  !is.na(hculture$V0[i]) & !is.na(hculture$V1[i]) & hculture$V0[i]==0 & hculture$V1[i]==0)
              {hculture$d_index[i]<-1}      else {hculture$d_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]=="A" | hculture$V1[i]=="B") ))
              {hculture$exd_contact[i]=1}   else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]   # for calculating secondary cases
d_index <- hculture[hculture$member==0,c(1,9)]     # for calculating SAR

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
hculture <- merge(hculture[,-8], d_index)
dim(hculture)

for (i in 1:nrow(hculture)){
    if ( hculture$exd_index[i]==1 | hculture$exd_contact[i] ==1)
       {hculture$exclude[i] <-1}
       else  {hculture$exclude[i] <-0}
}

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exclude[i] == 0 & !is.na(hculture$V1[i]) & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) | (hculture$V4[i] !=0 & !is.na(hculture$V4[i])) ) )
              {hculture$labsedcase[i] <- 1}
     else {hculture$labsedcase[i] <- 0}
}

###

visitday <- read.csv(paste(dir, "hchar_h.csv", sep=""), header=TRUE)
cdat <- read.csv(paste(dir, "clinicdat_h.csv", sep=""), header=TRUE)
lab <- hculture

visitday <- merge(visitday,cdat[,c(2,25)],by="hhID",all.x=TRUE)   # hhID, onsettime
visitday$onsettime <- floor(visitday$onsettime/2)
names(visitday)[9] <- "clinic_day"

lab2nd <- data.frame(hhID=unique(lab$hhID[lab$labsedcase==1]))
sedcase <- lab[lab$labsedcase==1,c(1,2)]

onset <- merge(visitday[,c(1,9)],lab2nd,by="hhID",all.y=TRUE)

#### -------------------------------------------------------------------------------------------------------------------------------------

symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""), header=TRUE)
symptom$day <- rep(0:9,nrow(symptom)/10)
symptom_2nd <- symptom[symptom$hhID==sedcase$hhID[1],]
for (i in 2:nrow(sedcase)){
     tmp <- symptom[symptom$hhID==sedcase$hhID[i],]
     symptom_2nd <- rbind(symptom_2nd,tmp)
}
symptom_2nd <- merge(symptom_2nd,onset,by="hhID",all.x=TRUE)         
symptom_2nd$day <- symptom_2nd$day + symptom_2nd$clinic_day

# set day 0 as index symtom onset, to obtain day of contact symtom onset
symptom_2nd$fever <- 1*(symptom_2nd$bodytemp >=37.8)
for (i in 1:nrow(symptom_2nd)){
      symptom_2nd$resp_symp[i] <- sum(symptom_2nd$sthroat[i],symptom_2nd$cough[i],
                                      symptom_2nd$pmuscle[i],symptom_2nd$headache[i],symptom_2nd$fever[i],na.rm=TRUE)
}


#
# Plot for h07075 m0&m2, h07117 m0&m4, h07132 m0&m4
#

p1 <- symptom_2nd[symptom_2nd$hhID=="h07075",]
p1_index <- p1[p1$member==0,]
p1_2nd <- p1[p1$member==2,]
p2 <- symptom_2nd[symptom_2nd$hhID=="h07117",]
p2_index <- p2[p2$member==0,]
p2_2nd <- p2[p2$member==4,]
p3 <- symptom_2nd[symptom_2nd$hhID=="h07132",]
p3 <- p3[1:50,]
p3_index <- p3[p3$member==0,]
p3_2nd <- p3[p3$member==4,]

windows(width=5.5, height=8)
layout(matrix(1:3, ncol=1))

# h07075 m2
par(mar=c(5,10,1.5,10), xaxs="i", yaxs="i")
plot(0, type="n", axes=FALSE, xlim=c(-0.5,10.5), ylim=c(0,11.8), xlab="", ylab="")

polygon(c(rep(c(0:10),each=2)), c(0,rep(c(p1_index$resp_symp),each=2),0), col=gray(0.7),lty=0)
axis(1, pos=0, cex.axis=1.0,at=0:10+0.5,labels=0:10)
axis(2, at=0:5,las=1,labels=0:5)
mtext("Day from index symptom onset", side=1,line=2.5,cex=0.8)
mtext("Number of
symptoms
(index case)", side=2,line=3,cex=0.8,at=2.5,las=1)

polygon(c(rep(c(0:10),each=2)), c(5.2,rep(p1_2nd$resp_symp+5.2,each=2),5.2), col=gray(0.7),lty=0)
axis(4, at=0:5+5.2,las=1,labels=0:5)
mtext("Number of 
symptoms
(secondary 
case)", side=4,line=3,cex=0.8,at=7.7,las=1)

lines(c(0,11.8),c(5.2,5.2))
lines(c(0,11.8),c(0,0))
lines(type= "p", c(1,3,7)+0.5, rep(10.8,3),lty=1,cex=1.2, pch=c("-","-","+"))
lines(c(3,3,7,7)+0.5,c(10.2,10.5,10.5,10.2))
title("(a)")
mtext("Sore
throat", side=1,line=-10.5,cex=0.6,at=-0.3)
arrows(0,8,0.5,6.5,length=0.05)

# h07117 m4
par(mar=c(5,10,1.5,10), xaxs="i", yaxs="i")
plot(0, type="n", axes=FALSE, xlim=c(-0.5,10.5), ylim=c(0,11.8), xlab="", ylab="")

polygon(c(rep(c(0:10),each=2)), c(0,rep(p2_index$resp_symp,each=2),0), col=gray(0.7),lty=0)
axis(1, pos=0, cex.axis=1.0,at=0:10+0.5,labels=0:10)
axis(2, at=0:5,las=1,labels=0:5)
mtext("Day from index symptom onset", side=1,line=2.5,cex=0.8)
mtext("Number of
symptoms
(index case)", side=2,line=3,cex=0.8,at=2.5,las=1)

polygon(c(rep(c(0:10),each=2)), c(5.2,rep(p2_2nd$resp_symp+5.2,each=2),5.2), col=gray(0.7),lty=0)
axis(4, at=0:5+5.2,las=1,labels=0:5)
mtext("Number of 
symptoms
(secondary 
case)", side=4,line=3,cex=0.8,at=7.7,las=1)

lines(c(0,11.8),c(5.2,5.2))
lines(c(0,11.8),c(0,0))
lines(type= "p", c(0,4,7,11)+0.5, rep(10.8,4),lty=1,cex=1.2, pch=c("-","-","+",NA))
lines(c(4,4,7,7)+0.5,c(10.2,10.5,10.5,10.2))
title("(b)")
mtext("Headache", side=1,line=-10.5,cex=0.6,at=-0.3)
arrows(0,8,0.5,6.5,length=0.05)

# h07132 m4
par(mar=c(5,10,1.5,10), xaxs="i", yaxs="i")
plot(0, type="n", axes=FALSE, xlim=c(-0.5,10.5), ylim=c(0,11.8), xlab="", ylab="")

polygon(c(rep(c(0:11),each=2)), c(0,rep(c(2,p3_index$resp_symp),each=2),0), col=gray(0.7),lty=0)
axis(1, pos=0, cex.axis=1.0,at=0:10+0.5,labels=0:10)
axis(2, at=0:5,las=1,labels=0:5)
mtext("Day from index symptom onset", side=1,line=2.5,cex=0.8)
mtext("Number of
symptoms
(index case)", side=2,line=3,cex=0.8,at=2.5,las=1)

polygon(c(rep(c(1:11),each=2)), c(5.2,rep(p3_2nd$resp_symp+5.2,each=2),5.2), col=gray(0.7),lty=0)
axis(4, at=0:5+5.2,las=1,labels=0:5)
mtext("Number of 
symptoms 
(secondary 
case)", side=4,line=3,cex=0.8,at=7.7,las=1)

lines(c(0,11.8),c(5.2,5.2))
lines(c(0,11.8),c(0,0))
lines(type= "p", c(1,3,5,9)+0.5, rep(10.8,4),lty=1,cex=1.2, pch=c("-","-","-","+"))
lines(c(5,5,9,9)+0.5,c(10.2,10.5,10.5,10.2))
title("(c)")
mtext("Sore
throat", side=1,line=-10.5,cex=0.6,at=1.2)
arrows(1.5,8,2,6.5,length=0.05)

# End of script

