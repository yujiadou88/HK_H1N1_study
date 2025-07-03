#
# R syntax to reproduce information for Appendix Table 1 from:
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

exd_index <- hculture[hculture$member==0,c(1,8)]   
d_index <- hculture[hculture$member==0,c(1,9)]     

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
visitday$v1_day <- visitday$v1_day+visitday$clinic_day
visitday$v2_day <- visitday$v2_day+visitday$clinic_day
visitday$v3_day <- visitday$v3_day+visitday$clinic_day
visitday$v4_day <- visitday$v4_day+visitday$clinic_day

lab2nd <- data.frame(hhID=unique(lab$hhID[lab$labsedcase==1]))
sedcase <- lab[lab$labsedcase==1,c(1,2)]

onset <- merge(visitday[,c(1,9)],lab2nd,by="hhID",all.y=TRUE)

#### -------------------------------------------------------------------------------------------------------------------------------------

symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""), header=TRUE)
symptom$day <- rep(0:9,nrow(symptom)/10)
symptom_2nd <- symptom[symptom$hhID==sedcase$hhID[1]&symptom$member==sedcase$member[1],]
for (i in 2:nrow(sedcase)){
     tmp <- symptom[symptom$hhID==sedcase$hhID[i]&symptom$member==sedcase$member[i],]
     symptom_2nd <- rbind(symptom_2nd,tmp)
}
symptom_2nd <- merge(symptom_2nd,onset,by="hhID",all.x=TRUE)         
symptom_2nd$day <- symptom_2nd$day + symptom_2nd$clinic_day

symptom_2nd$fever <- 1*(symptom_2nd$bodytemp >=37.8)
for (i in 1:nrow(symptom_2nd)){
      symptom_2nd$resp_symp[i] <- sum(symptom_2nd$sthroat[i],symptom_2nd$cough[i],
                                      symptom_2nd$pmuscle[i],symptom_2nd$headache[i],symptom_2nd$fever[i],na.rm=TRUE)
}
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07075"&symptom_2nd$day==0] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07117"&symptom_2nd$day==0] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07132"&symptom_2nd$day==1] <- 0
symptom_2nd$resp_symp[symptom_2nd$hhID=="h07132"&symptom_2nd$day==2] <- 0

symptom_2nd$mark <- 0
for (j in 1:nrow(sedcase)){
     i <- (j-1)*10+1
     endc <- j*10
     while(i <= endc){
        if (symptom_2nd$resp_symp[i]>=1) {
    symptom_2nd$mark[i] <- 1 
    break
    }
        else  i <- i+1
     }
}

# construct the table
sedcase <- merge(sedcase,onset,by="hhID",all.x=TRUE)
sedcase <- merge(sedcase,symptom_2nd[symptom_2nd$mark==1,1:3],by=c("hhID","member"),all.x=T)
names(sedcase)[4] <- "day_any_symp"

# Construct the 4th column of the table (lab confirmation time)
hculture <- hculture[hculture$labsedcase==1,c(1,2,4:7)]
hculture2 <- merge(hculture,visitday[c(1,5:8)],by="hhID",all.x=TRUE)
for (i in 1:nrow(hculture2)){
   if (hculture2$V2[i]!=0){
       hculture2$lab_L[i] <- hculture2$v1_day[i]
       hculture2$lab_R[i] <- hculture2$v2_day[i]
   }
   else if (hculture2$V3[i]!=0){
       hculture2$lab_L[i] <- hculture2$v2_day[i]
       hculture2$lab_R[i] <- hculture2$v3_day[i]
   }
   else if (hculture2$V4[i]!=0){
       hculture2$lab_L[i] <- hculture2$v3_day[i]
       hculture2$lab_R[i] <- hculture2$v4_day[i]
   }
}

sedcase <- cbind(sedcase,hculture2[11:12])
appt1 <- sedcase[order(sedcase$clinic_day,sedcase$day_any_symp),]
appt1[1] <- 1:21
names(appt1)[1] <- "ID"
appt1 <- appt1[-2]
appt1

# End of script
