#
# R syntax to reproduce information for Table 2 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al. 
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

tab2 <- matrix(rep(NA,54),ncol=6)

demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

# swine flu versus seasonal flu A (H1/H3)
swine <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
season <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$qPCR>0&(hc$AH1==1|hc$AH3==1)])

# clinical symptoms
symp$headache[is.na(symp$headache)] <- symp$sthroat[is.na(symp$sthroat)] <- symp$cough[is.na(symp$cough)] <-
symp$pmuscle[is.na(symp$pmuscle)] <- symp$rnose[is.na(symp$rnose)] <- symp$phlegm[is.na(symp$phlegm)] <- 0
symp$fever <- 1*(!is.na(symp$bodytemp)&symp$bodytemp>=37.8)
symp$flu1 <- 1*(symp$fever+symp$cough+symp$headache+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm>=2)  # ARI
symp$flu2 <- 1*(symp$fever==1 & (symp$cough+symp$sthroat>=1))  # ILI

## Function cluster bootstrap CI - PCR
cbCI <- function(data,nboot){
   idlist <- unique(data$hhID)
   sed <- sum(data$labsedcase)
   est <- sed/dim(data)[1]  # estimated SAR
   sar <- rep(NA,nboot)
   for (i in 1:nboot){
	 idsample <- sort(sample(idlist,replace=TRUE))
	 newdata <- data[data$hhID==idsample[1],]
	 for(j in 2:length(idsample)){
            temp <- data[data$hhID==idsample[j],]
	    newdata <- rbind(newdata,temp)
	 }
  	sar[i] <- sum(newdata$labsedcase)/dim(newdata)[1]
   }
   round(c(sed,est,quantile(sar,c(0.025,0.975))),2)
}
## End of function

## Function cluster bootstrap CI  - clinical def
cbCI2 <- function(data,nboot){
   idlist <- unique(data$hhID)
   sed <- sum(data$clinicsedcase)
   est <- sed/dim(data)[1]  # estimated SAR
   sar <- rep(NA,nboot)
   for (i in 1:nboot){
	 idsample <- sort(sample(idlist,replace=TRUE))
	 newdata <- data[data$hhID==idsample[1],]
	 for(j in 2:length(idsample)){
            temp <- data[data$hhID==idsample[j],]
	    newdata <- rbind(newdata,temp)
	 }
	sar[i] <- sum(newdata$clinicsedcase)/dim(newdata)[1]
   }
   round(c(sed,est,quantile(sar,c(0.025,0.975))),2)
}
## End of function

for (j in 1:2){
  if (j==1) {
    hc2 <-hc[hc$hhID %in% swine,]
    hc2$flu <- 1*(hc2$swine>0)
  }  
  if (j==2) {
    hc2 <-hc[hc$hhID %in% season,]
    hc2$flu <- 1*(hc2$qPCR>0)
  }
  
  hculture <- reshape(hc2[c(1:3,10)], timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="flu")
  names(hculture) <- c("hhID","member","V1","V2","V3")
  
  for (i in 1:nrow(hculture)){
  ## index exclusion: V1+V2+V3 culture is 0/NA
       if(hculture$member[i]==0 & !is.na(hculture$V1[i]) & hculture$V1[i]==0) {hculture$d_index[i]<-1}
         else {hculture$d_index[i]<-0}
       if(hculture$member[i]==0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]==0)|(is.na(hculture$V1[i])) )) {hculture$exd_index[i]<-1}
         else {hculture$exd_index[i]<-0}
       if(hculture$member[i]!=0 &  !is.na(hculture$V1[i]) & hculture$V1[i]!=0)
         {hculture$d_contact[i]<-1}
         else {hculture$d_contact[i]<-0}
       if(hculture$member[i]!=0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]!=0) | is.na(hculture$V1[i]) ))
            {hculture$exd_contact[i]=1}
  	  else{hculture$exd_contact[i]=0}
  }
  
  d_contactid <- unique(hculture$hhID[hculture$d_contact==1]) # for excluding co-index households
  d_contact <- data.frame(hhID=d_contactid)
  d_contact$d_contact <- 1
  
  exd_index <- hculture[hculture$member==0,c(1,7)]   # for excluding households with index -ve at baseline
  d_index <- hculture[hculture$member==0,c(1,6)]     # not used here
  
  #dim(hculture)
  hculture <- merge(hculture[,-7], exd_index)
  hculture <- merge(hculture[,-6], d_index)
  hculture <- merge(hculture[,-6], d_contact,all.x=TRUE)
  hculture$d_contact[is.na(hculture$d_contact)] <- 0
  #dim(hculture)
  hculture <- hculture[order(hculture$hhID,hculture$member),]
  
  hculture$analyzed <- 1*(hculture$exd_index==0&hculture$d_contact==0)
  
  # lab-2nd  
  for (i in 1:nrow(hculture)){
      if ( hculture$member[i] != 0 & hculture$analyzed[i]==1 & hculture$exd_contact[i]==0 &
               ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) | (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) ) )
           hculture$labsedcase[i] <- 1
  	  else hculture$labsedcase[i] <- 0
  }
  
  hculture.cc <- merge(hculture,hchar[c(1,3)],by="hhID",all.x=TRUE)
  names(hculture.cc)[dim(hculture.cc)[2]]<- "hhsize"
  hculture.cc2 <- hculture.cc[hculture.cc$member!=0&hculture.cc$analyzed==1&!is.na(hculture.cc$hhsize),]
  hculture.cc2$hhsize <- hculture.cc2$hhsize -1
  hculture.cc2$hhsize[hculture.cc2$member>1] <- 0
  hc_contact <- hculture.cc2[c(1,2,11,12)]
  
  set.seed(12345)
  tab2[1,1:3+(j-1)*3] <- cbCI(hc_contact,1000)[-1]
    
  # separate by index age group
  hc_contact2 <- merge(hc_contact,demog[demog$member==0,c(1,3)],by="hhID",all.x=TRUE)
  hc_contact_ic <- hc_contact2[!is.na(hc_contact2$age)&hc_contact2$age<=15,]
  hc_contact_ia <- hc_contact2[!is.na(hc_contact2$age)&hc_contact2$age>15,]
  set.seed(123456)
  tab2[4,1:3+(j-1)*3] <- cbCI(hc_contact_ic,1000)[-1]
  tab2[7,1:3+(j-1)*3] <- cbCI(hc_contact_ia,1000)[-1]
  
  
  
  # clinical SAR
  for(k in 1:2){
    hclinic <- data.frame(hhID = hculture$hhID)
    hclinic$member <- hculture$member
    symp.temp <- reshape(symp[c(1:3,11+k)], timevar="day", idvar=c("hhID","member"), direction="wide")
    hclinic <- merge(hclinic,symp.temp, by=c("hhID","member"), all.x=TRUE)
    hclinic <- hclinic[order(hclinic$hhID,hclinic$member),]
    names(hclinic) <- c("hhID","member","day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")
  
    hclinic$exd_contact <- hculture$exd_contact
    hclinic$analyzed <- hculture$analyzed
    
    ## Define secondary cases
    for (i in 1:nrow(hclinic)){
        if ( hclinic$member[i] != 0 & hclinic$analyzed[i] == 1 & hclinic$exd_contact[i]==0 &
             ( (!is.na(hclinic$day1[i]) & hclinic$day1[i]==1) | (!is.na(hclinic$day2[i]) & hclinic$day2[i]==1) |
    	   (!is.na(hclinic$day3[i]) & hclinic$day3[i]==1) | (!is.na(hclinic$day4[i]) & hclinic$day4[i]==1) |
    	   (!is.na(hclinic$day5[i]) & hclinic$day5[i]==1) | (!is.na(hclinic$day6[i]) & hclinic$day6[i]==1) |
    	   (!is.na(hclinic$day7[i]) & hclinic$day7[i]==1) | (!is.na(hclinic$day8[i]) & hclinic$day8[i]==1) |
    	   (!is.na(hclinic$day9[i]) & hclinic$day9[i]==1) ))
                  {hclinic$clinicsedcase[i] <- 1}
    	 else {hclinic$clinicsedcase[i] <- 0}
    }
    
    hclinic$hhsize <- hculture.cc$hhsize
    clicon <- hclinic[hclinic$member>0 & hclinic$analyzed==1 & !is.na(hclinic$hhsize),]
    clicon$hhsize <- hculture.cc2$hhsize
  
    set.seed(12345)
    tab2[1+k,1:3+(j-1)*3] <- cbCI2(clicon,1000)[-1]
    
    # separate by index age group
    clicon2 <- merge(clicon,demog[demog$member==0,c(1,3)],by="hhID",all.x=TRUE)
    clicon_ic <- clicon2[!is.na(clicon2$age)&clicon2$age<=15,]
    clicon_ia <- clicon2[!is.na(clicon2$age)&clicon2$age>15,]
    set.seed(123456)
    tab2[4+k,1:3+(j-1)*3] <- cbCI2(clicon_ic,1000)[-1]
    tab2[7+k,1:3+(j-1)*3] <- cbCI2(clicon_ia,1000)[-1]

  }
}

rownames(tab2) <- c("Any - RT-PCR","ARI","ILI","<16yr - RT-PCR","ARI","ILI",">16yr - RT-PCR","ARI","ILI")
colnames(tab2) <- c("pH1N1-SAR","95% CI - low","95% CI - up","sFluA-SAR","95% CI - low","95% CI - up")
tab2

# End of script