#
# R syntax to reproduce information for text results (Result section) from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al.
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

# Paragraph 1: ---------------------------------------------------------------------------------------------------------------------------------------
# QuickVue test vs gold standard (PCR)

clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

ph1 <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
sh1 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH1==1])
sh3 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH3==1])

ph1.all <- unique(c(as.character(clinic$scrID[clinic$hhID%in%ph1]),as.character(clinic$scrID[!is.na(clinic$pcr.pH1N1)&clinic$pcr.pH1N1==1])))
sh1.all <- unique(c(as.character(clinic$scrID[clinic$hhID%in%sh1]),as.character(clinic$scrID[!is.na(clinic$pcr.H1)&clinic$pcr.H1==1])))
sh3.all <- unique(c(as.character(clinic$scrID[clinic$hhID%in%sh3]),as.character(clinic$scrID[!is.na(clinic$pcr.H3)&clinic$pcr.H3==1])))
season.all <- c(sh1.all,sh3.all)
# function to calculate sensitivity
sens <- function(a,b){
    # a is true positive, b=total PCR+
    sens <- a/b
    CI <- binom.test(a,b)$conf[1:2]
    round(c(sens,CI),2)
}
# sensitivity (95% CI) for seasonal flu
sens(nrow(clinic[clinic$scrID%in%season.all&clinic$QVres==1,]),nrow(clinic[clinic$scrID%in%season.all,]))
# sensitivity (95% CI) for swine flu
sens(nrow(clinic[clinic$scrID%in%ph1.all&clinic$QVres==1,]),nrow(clinic[clinic$scrID%in%ph1.all,]))

# Paragraph 4: ---------------------------------------------------------------------------------------------------------------------------------------
# Serial intervals for pH1N1 and sH3N2

source("../Serial_scripts/functions.r")
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
ph1 <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
sh3 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH3==1])

# Define secondary cases first (by PCR confirmed)
hculture <- reshape(hc[c(1:4)], timevar="visit", idvar=c("hhID","member"), direction="wide")
names(hculture) <- c("hhID","member","V1","V2","V3")
for (i in 1:nrow(hculture)){
  ## index exclusion: V1+V2+V3 culture is 0/NA
     if(hculture$member[i]==0 & !is.na(hculture$V1[i]) & hculture$V1[i]==0) {hculture$d_index[i]<-1}
       else {hculture$d_index[i]<-0}
     if(hculture$member[i]==0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]==0)|(is.na(hculture$V1[i])) )) {hculture$exd_index[i]<-1}
       else {hculture$exd_index[i]<-0}
     if(hculture$member[i]!=0 &  !is.na(hculture$V1[i]) & hculture$V1[i]!=0)  {hculture$d_contact[i]<-1}
       else {hculture$d_contact[i]<-0}
     if(hculture$member[i]!=0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]!=0) | is.na(hculture$V1[i]) )) {hculture$exd_contact[i]=1}
	     else{hculture$exd_contact[i]=0}
}

d_contactid <- unique(hculture$hhID[hculture$d_contact==1]) # for excluding co-index households
d_contact <- data.frame(hhID=d_contactid)
d_contact$d_contact <- 1

exd_index <- hculture[hculture$member==0,c(1,7)]   # for calculating secondary cases
d_index <- hculture[hculture$member==0,c(1,6)]     # for excluding households with index -ve at baseline

hculture <- merge(hculture[,-7], exd_index)
hculture <- merge(hculture[,-6], d_index)
hculture <- merge(hculture[,-6], d_contact,all.x=TRUE)
hculture$d_contact[is.na(hculture$d_contact)] <- 0
hculture <- hculture[order(hculture$hhID,hculture$member),]
hculture$analyzed <- 1*(hculture$exd_index==0&hculture$d_contact==0)

# lab-2nd
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$analyzed[i]==1 & hculture$exd_contact[i]==0 &
             ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) | (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) ) )
         hculture$labsedcase[i] <- 1
	  else hculture$labsedcase[i] <- 0
}
sedcase <- hculture[hculture$labsedcase==1,c(1,2,11)]

# symtoms
symp$headache[is.na(symp$headache)] <- symp$sthroat[is.na(symp$sthroat)] <- symp$cough[is.na(symp$cough)] <-
symp$pmuscle[is.na(symp$pmuscle)] <- symp$rnose[is.na(symp$rnose)] <- symp$phlegm[is.na(symp$phlegm)] <- 0
symp$fever <- 1*(!is.na(symp$bodytemp)&symp$bodytemp>=37.8)
symp$flu <- 1*(symp$fever+symp$cough+symp$headache+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm>=2)

# get secondary cases only
symp2 <- merge(symp,sedcase,by=c("hhID","member"),all.x=TRUE)
symp2 <- symp2[symp2$labsedcase==1&!is.na(symp2$labsedcase),]
symp2 <- symp2[order(symp2$hhID,symp2$member,symp2$day),]

symp2$mark <- 0
for (j in 1:nrow(sedcase)){
     i <- (j-1)*10+1
     endc <- j*10
     while(i <= endc){
        if (symp2$flu[i]>=1) {
    symp2$mark[i] <- 1
    break
    }
        else  i <- i+1
     }
}

# set day 0 as index symptom onset day
sedcase <- merge(sedcase[1:2],symp2[symp2$mark==1,c("hhID","member","day")],by=c("hhID","member"),all.x=TRUE)
sedcase <- merge(sedcase,hchar[c(1,6,7)],all.x=TRUE)
sedcase$day[!is.na(sedcase$day)&sedcase$day==0] <- 1
sedcase$day <- sedcase$day+sedcase$v1_day

# Fit the weibull6 model (allow for truncted serial interval)

# for pH1N1 & sH3N2
sed.ph1 <- sedcase[sedcase$hhID%in%ph1,]
sed.sh3 <- sedcase[sedcase$hhID%in%sh3,]

pboot.weibull <- function(reps, time.b, data){
   # function to draw iid samples from fitted model and fit weibull models to each resample
  n <- length(data[,1])
  temp <- rep(NA, reps)
  output <- data.frame(mean=temp, par1=temp, par2=temp)
  fitted.weibull <- optim(c(0.5,1.5), weibull.loglik, time=time.b,data=data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
  for(i in 1:reps){
    new.data <- data.frame(rep(NA,2*n))
    new.data[,1] <- rweibull(2*n, exp(fitted.weibull$par[1]), exp(fitted.weibull$par[2]))
    names(new.data)[1] <- "day_symp"
    new.data$clinic_day <- sample(c(0,1,2),2*n,replace=TRUE,prob=c(7,13,1))             # proportion based on sample
    new.data <- new.data[new.data$day>=new.data$clinic_day,]
    new.data <- new.data[1:n,]
    new.weibull <- optim(c(0.5,1.5), weibull.loglik, time=new.data$day,data=new.data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
    output$mean[i] <- exp(new.weibull$par[2])*gamma(1+exp(-new.weibull$par[1]))
    output$par1[i] <- new.weibull$par[1]
    output$par2[i] <- new.weibull$par[2]
  }
  output
}

serial.weibull.ph1 <- optim(c(0.5,0.5), weibull.loglik, time=sed.ph1$day,data=sed.ph1, method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
mean.ph1 <- exp(serial.weibull.ph1$par[2])*gamma(1+exp(-serial.weibull.ph1$par[1]))
var.ph1 <- exp(2*serial.weibull.ph1$par[2])*(gamma(1+2*exp(-serial.weibull.ph1$par[1]))-(gamma(1+exp(-serial.weibull.ph1$par[1])))^2)
round(c(mean.ph1,sqrt(var.ph1)),1) # mean, sd
set.seed(12406)
wei.boot.ph1 <- pboot.weibull(500,sed.ph1$day,sed.ph1)
round(quantile(wei.boot.ph1$mean, c(0.025, 0.975)),1)    # 95% CI

serial.weibull.sh3 <- optim(c(0.5,1.5), weibull.loglik, time=sed.sh3$day,data=sed.sh3, method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
mean.sh3 <- exp(serial.weibull.sh3$par[2])*gamma(1+exp(-serial.weibull.sh3$par[1]))
var.sh3 <- exp(2*serial.weibull.sh3$par[2])*(gamma(1+2*exp(-serial.weibull.sh3$par[1]))-(gamma(1+exp(-serial.weibull.sh3$par[1])))^2)
round(c(mean.sh3,sqrt(var.sh3)),1) # mean, sd
set.seed(12406)
wei.boot.sh3 <- pboot.weibull(500,sed.sh3$day,sed.sh3)
round(quantile(wei.boot.sh3$mean, c(0.025, 0.975)),1)    # 95% CI

# Paragraph 9: ---------------------------------------------------------------------------------------------------------------------------------------
# Tamiflu vs antibody rise

sero <- read.csv(paste(dir, "paired_sera.csv", sep=""))
demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),1:9]

ph1 <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
swine <- unique(hc$hhID[hc$member==0&hc$visit==1&(hc$swine>0|hc$hasw>0)])
sh1 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH1==1])
sh3 <- unique(hc$hhID[hc$member==0&hc$visit==1&hc$AH3==1])

# index with H1N1 pandemic
pandemic <- sero[sero$hhID%in%swine&sero$member==0,]
pandemic$ratio.neut <- pandemic$neut_pH1_c/pandemic$neut_pH1
pandemic$ratio.hi <- pandemic$HI_pH1_c/pandemic$HI_pH1
pandemic <- merge(pandemic,demog[c(1,2,5)],by=c("hhID","member"),all.x=TRUE)

# HAI antibody titer rise (tamiflu vs no tamiflu)
round(c(exp(mean(log(pandemic$ratio.hi[pandemic$tamiflu==1]))),exp(mean(log(pandemic$ratio.hi[pandemic$tamiflu==0])))),1)
# p-value
round(wilcox.test(pandemic$ratio.hi[pandemic$tamiflu==1],pandemic$ratio.hi[pandemic$tamiflu==0])$p.value,2)
# p-value for difference in titer rise by viral neutralization
round(wilcox.test(pandemic$ratio.neut[pandemic$tamiflu==1],pandemic$ratio.neut[pandemic$tamiflu==0])$p.value,2)

# End of script
