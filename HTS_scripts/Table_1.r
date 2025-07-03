#
# R syntax to reproduce information for Table 1 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al. 
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

tab1 <- matrix(rep(NA,130),ncol=5)

clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

# swine flu versus seasonal flu A (H1/H3)
swine <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
season <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$qPCR>0&((hc$AH1==1&!is.na(hc$AH1))|(hc$AH3==1&!is.na(hc$AH3)))])

demog$agegp <- cut(demog$age,c(0,5,15,30,50,100))
demog.swi <- demog[demog$hhID%in%swine,]       # 45 hhs
demog.sea <- demog[demog$hhID%in%season,]    # 54 hhs
clinic.swi <- clinic[clinic$hhID%in%swine,]
clinic.sea <- clinic[clinic$hhID%in%season,]

# index characteristics
index <- c(nrow(demog.swi[demog.swi$member==0,]),nrow(demog.sea[demog.sea$member==0,]))
index
tab1[1:5,1] <- table(demog.swi$agegp[demog.swi$member==0])
tab1[6,1] <- sum(demog.swi$male[demog.swi$member==0]) 
tab1[7:13,1] <- colSums(clinic.swi[c(6:12)])
tab1[14,1] <- sum(demog.swi$tamiflu[demog.swi$member==0])
tab1[15:18,1] <- table(clinic.swi$onsettime)
tab1[19,1] <- sum(demog.swi$vaccine09[demog.swi$member==0])
tab1[1:5,3] <- table(demog.sea$agegp[demog.sea$member==0])
tab1[6,3] <- sum(demog.sea$male[demog.sea$member==0])
tab1[7:13,3] <- colSums(clinic.sea[c(6:12)],na.rm=T)
tab1[14,3] <- sum(demog.sea$tamiflu[demog.sea$member==0])
tab1[15:18,3] <- table(clinic.sea$onsettime)
tab1[19,3] <- sum(demog.sea$vaccine09[demog.sea$member==0])
tab1[1:19,2] <- round(tab1[1:19,1]/index[1],2)
tab1[1:19,4] <- round(tab1[1:19,3]/index[2],2)
tab1[5,5] <- round(fisher.test(tab1[1:5,c(1,3)])$p.value,2)
for (i in 6:14){
  tab1[i,5] <- round(chisq.test(rbind(tab1[i,c(1,3)],index-tab1[i,c(1,3)]))$p.value,2)
}
tab1[18,5] <- round(chisq.test(tab1[15:18,c(1,3)])$p.value,2)
tab1[19,5] <- round(chisq.test(rbind(tab1[19,c(1,3)],index-tab1[19,c(1,3)]))$p.value,2)

# contact characteristics
contact <- c(nrow(demog.swi[demog.swi$member>0,]),nrow(demog.sea[demog.sea$member>0,]))
contact
tab1[20:24,1] <- table(demog.swi$agegp[demog.swi$member>0])
tab1[25,1] <- sum(demog.swi$male[demog.swi$member>0])
tab1[26,1] <- sum(demog.swi$vaccine09[demog.swi$member>0],na.rm=T)
tab1[20:24,3] <- table(demog.sea$agegp[demog.sea$member>0])
tab1[25,3] <- sum(demog.sea$male[demog.sea$member>0])
tab1[26,3] <- sum(demog.sea$vaccine09[demog.sea$member>0])
tab1[20:26,2] <- round(tab1[20:26,1]/contact[1],2)
tab1[20:26,4] <- round(tab1[20:26,3]/contact[2],2)
tab1[24,5] <- round(fisher.test(tab1[20:24,c(1,3)])$p.value,2)
tab1[25,5] <- round(chisq.test(rbind(tab1[25,c(1,3)],contact-tab1[25,c(1,3)]))$p.value,2)
tab1[26,5] <- round(chisq.test(rbind(tab1[26,c(1,3)],contact-tab1[26,c(1,3)]))$p.value,2)

colnames(tab1) <- c("pH1N1","%","Seasonal flu A","%","p-value")
rownames(tab1) <- c("Index age <=5yr","6-15yr","16-30yr","31-50yr",">50yr","Male","Fever","Headache","Sore throat","Cough","Myalgia","Runny nose",
                    "Phlegm","Tamiflu","0-12hr","13-24hr","25-36hr","37-48hr","Vax08-09",
                    "Contact age <=5yr","6-15yr","16-30yr","31-50yr",">50yr","Male","Vax08-09")
tab1

# End of script

