#
# R syntax to reproduce information for Appendix Table 2 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al.
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

appt2 <- matrix(rep(NA,80),ncol=8)

symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
sero <- read.csv(paste(dir, "paired_sera.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

swine <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
sh3 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH3==1&!is.na(hc$AH3)])

sero$rise.ph1 <- 1*(sero$neut_pH1_c/sero$neut_pH1>=4)
sero$rise.sh3 <- 1*(sero$neut_sH3_c/sero$neut_sH3>=4)

pcr <- unique(hc[hc$member>0,1:2])

# pH1N1 - PCR
contact.ph1 <- sero[sero$hhID%in%swine&sero$member>0,]
pcr.ph1 <- unique(hc[hc$swine==1&!is.na(hc$swine),1:2])
pcr.ph1$pcr.ph1 <- 1
pcr2 <- merge(pcr,pcr.ph1,by=c("hhID","member"),all.x=TRUE)
pcr2$pcr.ph1[is.na(pcr2$pcr.ph1)] <- 0
contact.ph1 <- merge(contact.ph1,pcr2,by=c("hhID","member"),all.x=T)
contact.ph1 <- contact.ph1[!is.na(contact.ph1$pcr.ph1)&!is.na(contact.ph1$rise.ph1),]

appt2[1,1] <- nrow(contact.ph1[contact.ph1$rise.ph1==0&contact.ph1$pcr.ph1==0,])
appt2[1,3] <- nrow(contact.ph1[contact.ph1$rise.ph1==0&contact.ph1$pcr.ph1==1,])
appt2[1,c(5,7)] <- table(contact.ph1$pcr.ph1[contact.ph1$rise.ph1==1])

# sH3N2 - PCR
contact.sh3 <- sero[sero$hhID%in%sh3&sero$member>0,]
pcr.sh3 <- unique(hc[hc$AH3==1&!is.na(hc$AH3),1:2])
pcr.sh3$pcr.sh3 <- 1
pcr3 <- merge(pcr,pcr.sh3,by=c("hhID","member"),all.x=TRUE)
pcr3$pcr.sh3[is.na(pcr3$pcr.sh3)] <- 0
contact.sh3 <- merge(contact.sh3,pcr3,by=c("hhID","member"),all.x=T)
contact.sh3 <- contact.sh3[!is.na(contact.sh3$pcr.sh3)&!is.na(contact.sh3$rise.sh3),]

appt2[6,1] <- nrow(contact.sh3[contact.sh3$rise.sh3==0&contact.sh3$pcr.sh3==0,])
appt2[6,3] <- nrow(contact.sh3[contact.sh3$rise.sh3==0&contact.sh3$pcr.sh3==1,])
appt2[6,c(5,7)] <- table(contact.sh3$pcr.sh3[contact.sh3$rise.sh3==1])

## symptoms
symp$fever <- 1*(symp$bodytemp>=37.8)
symp$fever[is.na(symp$fever)] <- symp$headache[is.na(symp$headache)] <- symp$sthroat[is.na(symp$sthroat)] <-
symp$cough[is.na(symp$cough)] <- symp$pmuscle[is.na(symp$pmuscle)] <- symp$rnose[is.na(symp$rnose)] <-
symp$phlegm[is.na(symp$phlegm)] <- 0

symp$ARI <- 1*(symp$fever+symp$cough+symp$headache+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm>=2)
symp$ILI <- 1*(symp$fever==1 & symp$cough+symp$sthroat>=1)

sympf <- reshape(symp[c(1:3,11)],timevar="day",idvar=c("hhID","member"), direction="wide", v.names="fever")
for(i in 1:nrow(sympf)){
  sympf$fever[i] <- 1*(sum(sympf[i,3:12])>=1)
}
sympc <- reshape(symp[c(1:3,7)],timevar="day",idvar=c("hhID","member"), direction="wide", v.names="cough")
for(i in 1:nrow(sympc)){
  sympc$cough[i] <- 1*(sum(sympc[i,3:12])>=1)
}
sympARI <- reshape(symp[c(1:3,12)],timevar="day",idvar=c("hhID","member"), direction="wide", v.names="ARI")
for(i in 1:nrow(sympARI)){
  sympARI$ARI[i] <- 1*(sum(sympARI[i,3:12])>=1)
}
sympILI <- reshape(symp[c(1:3,13)],timevar="day",idvar=c("hhID","member"), direction="wide", v.names="ILI")
for(i in 1:nrow(sympILI)){
  sympILI$ILI[i] <- 1*(sum(sympILI[i,3:12])>=1)
}

# number of symptoms
rise.ph1 <- merge(contact.ph1[c(1:2,13,15)],sympf[c("hhID","member","fever")],by=c("hhID","member"),all.x=TRUE)
rise.ph1 <- merge(rise.ph1,sympc[c("hhID","member","cough")],by=c("hhID","member"),all.x=TRUE)
rise.ph1 <- merge(rise.ph1,sympARI[c("hhID","member","ARI")],by=c("hhID","member"),all.x=TRUE)
rise.ph1 <- merge(rise.ph1,sympILI[c("hhID","member","ILI")],by=c("hhID","member"),all.x=TRUE)

rise.sh3 <- merge(contact.sh3[c(1:2,14,15)],sympf[c("hhID","member","fever")],by=c("hhID","member"),all.x=TRUE)
rise.sh3 <- merge(rise.sh3,sympc[c("hhID","member","cough")],by=c("hhID","member"),all.x=TRUE)
rise.sh3 <- merge(rise.sh3,sympARI[c("hhID","member","ARI")],by=c("hhID","member"),all.x=TRUE)
rise.sh3 <- merge(rise.sh3,sympILI[c("hhID","member","ILI")],by=c("hhID","member"),all.x=TRUE)

appt2[2:5,1] <- colSums(rise.ph1[rise.ph1$rise.ph1==0&rise.ph1$pcr.ph1==0,5:8])
appt2[2:5,5] <- colSums(rise.ph1[rise.ph1$rise.ph1==1&rise.ph1$pcr.ph1==0,5:8])
appt2[2:5,7] <- colSums(rise.ph1[rise.ph1$rise.ph1==1&rise.ph1$pcr.ph1==1,5:8])
appt2[7:10,1] <- colSums(rise.sh3[rise.sh3$rise.sh3==0&rise.sh3$pcr.sh3==0,5:8])
appt2[7:10,5] <- colSums(rise.sh3[rise.sh3$rise.sh3==1&rise.sh3$pcr.sh3==0,5:8])
appt2[7:10,7] <- colSums(rise.sh3[rise.sh3$rise.sh3==1&rise.sh3$pcr.sh3==1,5:8])

appt2[2:5,2] <- round(appt2[2:5,1]/appt2[1,1],2)
appt2[2:5,6] <- round(appt2[2:5,5]/appt2[1,5],2)
appt2[2:5,8] <- round(appt2[2:5,7]/appt2[1,7],2)
appt2[7:10,2] <- round(appt2[7:10,1]/appt2[6,1],2)
appt2[7:10,6] <- round(appt2[7:10,5]/appt2[6,5],2)
appt2[7:10,8] <- round(appt2[7:10,7]/appt2[6,7],2)

colnames(appt2) <- c("ABrise-,PCR-","%","ABrise-,PCR+","%","ABrise+,PCR-","%","ABrise+,PCR+","%")
rownames(appt2) <- c("pH1N1","Fever","Cough","ARI","ILI","sH3N2","Fever","Cough","ARI","ILI")

appt2

# End of script
