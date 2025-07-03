#
# R syntax to reproduce information for Figure 2 from:
#
# Cowling BJ, Chan KH, Fang VJ, Lau LLH, So THC, et al.
# Comparative epidemiology of pandemic and seasonal influenza A in households
# NEJM, 2010 (in press).
#
# Last updated by Fang VJ, Lau LLH and Cowling BJ.
# April 7, 2010

dir <- "../data/HongKongHTSV1/"

geometric <- function(x)exp(mean(log(x)))
sero <- read.csv(paste(dir, "paired_sera.csv", sep=""))
demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
hchar <- read.csv(paste(dir, "hchar_h.csv", sep="")) 
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hc <- hc[!is.na(hc$qPCR),]

lab <- hc[c("hhID","member","visit","swine", "qPCR")]
lab <- reshape(lab, timevar="visit", idvar=c("hhID", "member"), direction="wide")
lab$pos <- lab$swine.1 + lab$swine.2 + lab$swine.3
lab$pos <- (lab$pos > 0)*1
lab$pos[is.na(lab$pos)] <- 0
lab$pos[lab$hhID == 751 & lab$member == 3] <- 1

# set fever
symp$fever <- 1*(symp$bodytemp >= 37.8)

###                                                                                 ###
###                                      SWINE FLU                                  ###
###                                                                                 ###

# isolate swine flu positives
swine <- lab[lab$pos == 1, ]

# set detection limit at 450
swine$qPCR.1[swine$qPCR.1 <= 900 & !is.na(swine$qPCR.1)] <- 450
swine$qPCR.2[swine$qPCR.2 <= 900 & !is.na(swine$qPCR.2)] <- 450
swine$qPCR.3[swine$qPCR.3 <= 900 & !is.na(swine$qPCR.3)] <- 450

# add demographic info
swine <- merge(swine, demog[1:4], by=c("hhID", "member"), all.x=TRUE)
swine$adult <- (swine$age >= 16)*1

# reorganize
swine <- swine[c(1:3,5,7,4,6,8,11,12)]

# isolate swine index subjects
swine_i <- swine[swine$member == 0, ]
swine_i <- merge(swine_i, hchar[c("hhID","v1_day","v2_day","v3_day")], by="hhID")
swine_i$v1_day <- as.numeric(as.character(swine_i$v1_day))
swine_i$v2_day <- as.numeric(as.character(swine_i$v2_day))
swine_i$v3_day <- as.numeric(as.character(swine_i$v3_day))

plot_i <- data.frame(hhID = rep(NA, (dim(swine_i)[1])*3))
plot_i$hhID <- c(swine_i$hhID, swine_i$hhID, swine_i$hhID)
plot_i$pcr <- c(swine_i$qPCR.1, swine_i$qPCR.2, swine_i$qPCR.3)
plot_i$day <- c(swine_i$v1_day, swine_i$v2_day, swine_i$v3_day)
plot_i$adult <- c(swine_i$adult, swine_i$adult, swine_i$adult)
plot_i <- plot_i[plot_i$pcr > 0, ]
plot_i$pcr <- log10(plot_i$pcr)

# mean shedding by day
mean_i <- data.frame(day = rep(NA, 10))
for(i in 1:10){
  mean_i$day[i] <- (i-1)
  geo <- plot_i$pcr[plot_i$day==(i-1)]
  mean_i$gmean[i] <- mean(geo)
}

# swine index symptoms
sw_symp_i <- swine_i[c(1,2)]
sw_symp_i <- merge(sw_symp_i, symp)
sw_symp_i <- sw_symp_i[order(sw_symp_i$hhID, sw_symp_i$member, sw_symp_i$day), ]

for (i in 1:dim(sw_symp_i)[1]){
    sw_symp_i$up[i] <- sum(sw_symp_i$sthroat[i], sw_symp_i$rnose[i], na.rm=TRUE)
    sw_symp_i$lo[i] <- sum(sw_symp_i$cough[i], sw_symp_i$phlegm[i], na.rm=TRUE)
    sw_symp_i$sy[i] <- sum(sw_symp_i$fever[i], sw_symp_i$headache[i], sw_symp_i$pmuscle[i], na.rm=TRUE)
    }

# readjust day by v1_day
sw_symp_i <- merge(sw_symp_i, swine_i[c(1,11)], by="hhID")
sw_symp_i$day <- sw_symp_i$day + sw_symp_i$v1_day

    
# swine index symptom plot
sw_i_symp_plot <- data.frame(day=rep(NA, 12))

for (i in 1:12) {
sw_i_symp_plot$day[i] <- (i-1)
sw_i_symp_plot$up[i] <- mean(sw_symp_i$up[sw_symp_i$day == (i-1)], na.rm=TRUE)
sw_i_symp_plot$lo[i] <- mean(sw_symp_i$lo[sw_symp_i$day == (i-1)], na.rm=TRUE)
sw_i_symp_plot$sy[i] <- mean(sw_symp_i$sy[sw_symp_i$day == (i-1)], na.rm=TRUE)
sw_i_symp_plot$temp[i] <- mean(sw_symp_i$bodytemp[sw_symp_i$day == (i-1)], na.rm=TRUE)
}

sw_i_symp_plot$up <- sw_i_symp_plot$up/2
sw_i_symp_plot$lo <- sw_i_symp_plot$lo/2
sw_i_symp_plot$sy <- sw_i_symp_plot$sy/3




# isolate secondary swine subjects
swine_s <- swine[swine$member != 0, ]
swine_s <- merge(swine_s, hchar[c("hhID","v1_day","v2_day","v3_day")], by="hhID")
swine_s$v1_day <- as.numeric(as.character(swine_s$v1_day))
swine_s$v2_day <- as.numeric(as.character(swine_s$v2_day))
swine_s$v3_day <- as.numeric(as.character(swine_s$v3_day))
swine_s <- swine_s[swine_s$hhID != 785,] # exclude swine +ve but 0 PCR

swine_s <- swine_s[swine_s$swine.1 == 0,] 


# finding ARI onset date
symp_s <- swine_s[c(1,2)]
symp_s <- merge(symp_s, symp, by=c("hhID","member"), all.x=TRUE)
symp_s <- symp_s[order(symp_s$hhID, symp_s$member, symp_s$day),]

for(i in 1:dim(symp_s)[1]){
symp_s$tot_symp[i] <- sum(symp_s$headache[i], symp_s$sthroat[i], symp_s$cough[i], symp_s$pmuscle[i], symp_s$rnose[i], 
                          symp_s$phlegm[i], symp_s$fever[i], na.rm=TRUE)
}

symp_s <- symp_s[c(1:3,12)]
symp_s <- reshape(symp_s, timevar="day", idvar=c("hhID", "member"), direction="wide")

symp_s$ari_onset <- 0 ### so that asymptomatics will be at day 0
symp_s$ari_onset[symp_s$tot_symp.0 >= 2] <- 0
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 >= 2] <- 1
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 <= 1 & symp_s$tot_symp.2 >= 2] <- 2
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 <= 1 & symp_s$tot_symp.2 <= 1 
          &symp_s$tot_symp.3 >= 2 ] <- 3
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 <= 1 & symp_s$tot_symp.2 <= 1 
          &symp_s$tot_symp.3 <= 1 & symp_s$symp.4 >= 2] <- 4
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 <= 1 & symp_s$tot_symp.2 <= 1 
          &symp_s$tot_symp.3 <= 1 & symp_s$symp.4 <= 1 & symp_s$sym.5 >= 2] <- 5
symp_s$ari_onset[symp_s$tot_symp.0 <= 1 & symp_s$tot_symp.1 <= 1 & symp_s$tot_symp.2 <= 1 
          &symp_s$tot_symp.3 <= 1 & symp_s$symp.4 <= 1 & symp_s$symp.5 <= 1 & symp_s$symp.6 >= 2] <- 6


# correct day accoring to ari onset
swine_s$v3_day <- swine_s$v3_day - symp_s$ari_onset - swine_s$v1_day
swine_s$v2_day <- swine_s$v2_day - symp_s$ari_onset - swine_s$v1_day
swine_s$v1_day <- swine_s$v1_day - symp_s$ari_onset - swine_s$v1_day

# remove asymptomatics, subclinical and no onset
swine_s <- swine_s[swine_s$v1_day != 0,]


plot_s <- data.frame(hhID = rep(NA, (dim(swine_s)[1])*3))
plot_s$hhID <- c(swine_s$hhID, swine_s$hhID, swine_s$hhID)
plot_s$pcr <- c(swine_s$qPCR.1, swine_s$qPCR.2, swine_s$qPCR.3)
plot_s$day <- c(swine_s$v1_day, swine_s$v2_day, swine_s$v3_day)
plot_s$adult <- c(swine_s$adult, swine_s$adult, swine_s$adult)
plot_s <- plot_s[plot_s$pcr > 0, ]
plot_s$pcr <- log10(plot_s$pcr)
plot_s <- plot_s[!is.na(plot_s$hhID),]

# mean shedding by day
mean_s <- data.frame(day = rep(NA, 14))
for(i in 1:14){
  mean_s$day[i] <- (i-4)
  geo <- plot_s$pcr[plot_s$day==(i-4)]
  mean_s$gmean[i] <- mean(geo)
}

# swine secondary symptoms
sw_symp_s <- swine_s[c(1,2)]
sw_symp_s <- merge(sw_symp_s, symp)
sw_symp_s <- sw_symp_s[order(sw_symp_s$hhID, sw_symp_s$member, sw_symp_s$day), ]

for (i in 1:dim(sw_symp_s)[1]){
    sw_symp_s$up[i] <- sum(sw_symp_s$sthroat[i], sw_symp_s$rnose[i], na.rm=TRUE)
    sw_symp_s$lo[i] <- sum(sw_symp_s$cough[i], sw_symp_s$phlegm[i], na.rm=TRUE)
    sw_symp_s$sy[i] <- sum(sw_symp_s$fever[i], sw_symp_s$headache[i], sw_symp_s$pmuscle[i], na.rm=TRUE)
    }
    
# add date of ARI onset
sw_symp_s <- merge(sw_symp_s, symp_s[c(1,2,13)], by=c("hhID","member"))
sw_symp_s$day <- sw_symp_s$day - sw_symp_s$ari_onset

# swine secondary symptom plot
sw_s_symp_plot <- data.frame(day=rep(NA, 14))

for (i in 1:14) {
sw_s_symp_plot$day[i] <- (i-4)
sw_s_symp_plot$up[i] <- mean(sw_symp_s$up[sw_symp_s$day == (i-4)], na.rm=TRUE)
sw_s_symp_plot$lo[i] <- mean(sw_symp_s$lo[sw_symp_s$day == (i-4)], na.rm=TRUE)
sw_s_symp_plot$sy[i] <- mean(sw_symp_s$sy[sw_symp_s$day == (i-4)], na.rm=TRUE)
sw_s_symp_plot$temp[i] <- mean(sw_symp_s$bodytemp[sw_symp_s$day == (i-4)], na.rm=TRUE)
}

sw_s_symp_plot$up <- sw_s_symp_plot$up/2
sw_s_symp_plot$lo <- sw_s_symp_plot$lo/2
sw_s_symp_plot$sy <- sw_s_symp_plot$sy/3

### TCID50 qculture #
tcid <- hc[c("hhID", "member", "visit", "q_culture")]
tcid <- tcid[order(tcid$hhID, tcid$q_culture), ]
tcid <- reshape(tcid, timevar="visit", idvar=c("hhID", "member"), direction="wide")
tcid$q_culture.1[tcid$q_culture.1 == 0] <- 0.3
tcid$q_culture.2[tcid$q_culture.2 == 0] <- 0.3
tcid$q_culture.3[tcid$q_culture.3 == 0] <- 0.3

# swine index tcid50
swine_i_tcid <- swine_i[c(1,2,10:13)]
swine_i_tcid <- merge(swine_i_tcid, tcid, by=c("hhID", "member"))

# swine index tcid50 plot
si_tcid_plot <- data.frame(day=rep(NA,10), tcid=rep(NA,10))

for (i in 1:10){
  si_tcid_plot[i,1] <- i
  si_tcid_plot[i,2] <- mean(c(swine_i_tcid$q_culture.1[swine_i_tcid$v1_day == (i)], swine_i_tcid$q_culture.2[swine_i_tcid$v2_day == (i)], swine_i_tcid$q_culture.3[swine_i_tcid$v3_day == (i)]), na.rm=TRUE)
}

# swine secondary tcid50
swine_s_tcid <- swine_s[c(1,2,10:13)]
swine_s_tcid <- merge(swine_s_tcid, tcid, by=c("hhID", "member"))

# swine sceondary tcid50 plot
ss_tcid_plot <- data.frame(day=rep(NA,11), tcid=rep(NA,11))

for (i in 1:11){
  ss_tcid_plot[i,1] <- (i-4)
  ss_tcid_plot[i,2] <- mean(c(swine_s_tcid$q_culture.1[swine_s_tcid$v1_day == (i-4)], swine_s_tcid$q_culture.2[swine_s_tcid$v2_day == (i-4)], swine_s_tcid$q_culture.3[swine_s_tcid$v3_day == (i-4)]), na.rm=TRUE)
}


###                                                                              ###
###                               SEASONAL FLU                                   ###
###                                                                              ###


seasonal <- lab[lab$pos != 1, ]
seasonal <- seasonal[c(1,2,4,6,8)]
seasonal$pos <- seasonal$qPCR.1 + seasonal$qPCR.2 + seasonal$qPCR.3
seasonal$pos <- (seasonal$pos > 0)*1
seasonal <- seasonal[seasonal$pos == 1 & !is.na(seasonal$pos),]

# set detection limit
seasonal$qPCR.1[seasonal$qPCR.1 <= 900 & !is.na(seasonal$qPCR.1)] <- 450
seasonal$qPCR.2[seasonal$qPCR.2 <= 900 & !is.na(seasonal$qPCR.2)] <- 450
seasonal$qPCR.3[seasonal$qPCR.3 <= 900 & !is.na(seasonal$qPCR.3)] <- 450

# add demographic info
seasonal <- merge(seasonal,  demog[1:4], by=c("hhID", "member"), all.x=TRUE)
seasonal$adult <- (seasonal$age >= 16)*1

# reorganize
seasonal <- seasonal[c(1:5,8:9)]

# isolate seasonal index subjects
seasonal_i <- seasonal[seasonal$member == 0, ]
seasonal_i <- merge(seasonal_i, hchar[c("hhID","v1_day","v2_day","v3_day")], by="hhID")
seasonal_i$v1_day <- as.numeric(as.character(seasonal_i$v1_day))
seasonal_i$v2_day <- as.numeric(as.character(seasonal_i$v2_day))
seasonal_i$v3_day <- as.numeric(as.character(seasonal_i$v3_day))

plot_seasonal_i <- data.frame(hhID = rep(NA, (dim(seasonal_i)[1])*3))
plot_seasonal_i$hhID <- c(seasonal_i$hhID, seasonal_i$hhID, seasonal_i$hhID)
plot_seasonal_i$pcr <- c(seasonal_i$qPCR.1, seasonal_i$qPCR.2, seasonal_i$qPCR.3)
plot_seasonal_i$day <- c(seasonal_i$v1_day, seasonal_i$v2_day, seasonal_i$v3_day)
plot_seasonal_i$adult <- c(seasonal_i$adult, seasonal_i$adult, seasonal_i$adult)
plot_seasonal_i <- plot_seasonal_i[plot_seasonal_i$pcr > 0, ]
plot_seasonal_i$pcr <- log10(plot_seasonal_i$pcr)

# mean shedding by day
mean_seasonal_i <- data.frame(day = rep(NA, 10))
for(i in 1:10){
  mean_seasonal_i$day[i] <- (i-1)
  geo <- plot_seasonal_i$pcr[plot_seasonal_i$day==(i-1)]
  mean_seasonal_i$gmean[i] <- mean(geo)
}

# seasonal index symptoms
se_symp_i <- seasonal_i[c(1,2)]
se_symp_i <- merge(se_symp_i, symp)
se_symp_i <- se_symp_i[order(se_symp_i$hhID, se_symp_i$member, se_symp_i$day), ]

for (i in 1:dim(se_symp_i)[1]){
    se_symp_i$up[i] <- sum(se_symp_i$sthroat[i], se_symp_i$rnose[i], na.rm=TRUE)
    se_symp_i$lo[i] <- sum(se_symp_i$cough[i], se_symp_i$phlegm[i], na.rm=TRUE)
    se_symp_i$sy[i] <- sum(se_symp_i$fever[i], se_symp_i$headache[i], se_symp_i$pmuscle[i], na.rm=TRUE)
}
  
# readjust day by v1_day
se_symp_i <- merge(se_symp_i, seasonal_i[c(1,8)], by="hhID")
se_symp_i$day <- se_symp_i$day + se_symp_i$v1_day
  
# seasonal index symptom plot
se_i_symp_plot <- data.frame(day=rep(NA, 12))

for (i in 1:12) {
se_i_symp_plot$day[i] <- (i-1)
se_i_symp_plot$up[i] <- mean(se_symp_i$up[se_symp_i$day == (i-1)], na.rm=TRUE)
se_i_symp_plot$lo[i] <- mean(se_symp_i$lo[se_symp_i$day == (i-1)], na.rm=TRUE)
se_i_symp_plot$sy[i] <- mean(se_symp_i$sy[se_symp_i$day == (i-1)], na.rm=TRUE)
se_i_symp_plot$temp[i] <- mean(se_symp_i$bodytemp[se_symp_i$day == (i-1)], na.rm=TRUE)
}

se_i_symp_plot$up <- se_i_symp_plot$up/2
se_i_symp_plot$lo <- se_i_symp_plot$lo/2
se_i_symp_plot$sy <- se_i_symp_plot$sy/3

# isolate seasonal secondary subjects
seasonal_s <- seasonal[seasonal$member != 0, ]
seasonal_s <- merge(seasonal_s, hchar[c("hhID","v1_day","v2_day","v3_day")], by="hhID")
seasonal_s$v1_day <- as.numeric(as.character(seasonal_s$v1_day))
seasonal_s$v2_day <- as.numeric(as.character(seasonal_s$v2_day))
seasonal_s$v3_day <- as.numeric(as.character(seasonal_s$v3_day))
seasonal_s <- seasonal_s[seasonal_s$qPCR.1 == 450,] 

# finding ARI onset date
symp_seasonal_s <- seasonal_s[c(1,2)]
symp_seasonal_s <- merge(symp_seasonal_s, symp, by=c("hhID","member"), all.x=TRUE)
symp_seasonal_s <- symp_seasonal_s[order(symp_seasonal_s$hhID, symp_seasonal_s$member, symp_seasonal_s$day),]

for(i in 1:dim(symp_seasonal_s)[1]){
symp_seasonal_s$tot_symp[i] <- sum(symp_seasonal_s$headache[i], symp_seasonal_s$sthroat[i], symp_seasonal_s$cough[i],
                                   symp_seasonal_s$pmuscle[i], symp_seasonal_s$rnose[i], symp_seasonal_s$phlegm[i], symp_seasonal_s$fever[i], na.rm=TRUE)
}

symp_seasonal_s <- symp_seasonal_s[c(1:3,12)]
symp_seasonal_s <- reshape(symp_seasonal_s, timevar="day", idvar=c("hhID", "member"), direction="wide")

symp_seasonal_s$ari_onset <- 0 ### so that asymptomatics will be at day 0
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 >= 2] <- 0
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 >= 2] <- 1
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 <= 1 & symp_seasonal_s$tot_symp.2 >= 2] <- 2
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 <= 1 & symp_seasonal_s$tot_symp.2 <= 1 
          &symp_seasonal_s$tot_symp.3 >= 2 ] <- 3
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 <= 1 & symp_seasonal_s$tot_symp.2 <= 1 
          &symp_seasonal_s$tot_symp.3 <= 1 & symp_seasonal_s$tot_symp.4 >= 2] <- 4
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 <= 1 & symp_seasonal_s$tot_symp.2 <= 1 
          &symp_seasonal_s$tot_symp.3 <= 1 & symp_seasonal_s$tot_symp.4 <= 1 & symp_seasonal_s$sym.5 >= 2] <- 5
symp_seasonal_s$ari_onset[symp_seasonal_s$tot_symp.0 <= 1 & symp_seasonal_s$tot_symp.1 <= 1 & symp_seasonal_s$tot_symp.2 <= 1 
          &symp_seasonal_s$tot_symp.3 <= 1 & symp_seasonal_s$tot_symp.4 <= 1 & symp_seasonal_s$symp.5 <= 1 & symp_seasonal_s$symp.6 >= 2] <- 6

# correct day accoring to ari onset
seasonal_s$v3_day <- seasonal_s$v3_day - symp_seasonal_s$ari_onset - seasonal_s$v1_day
seasonal_s$v2_day <- seasonal_s$v2_day - symp_seasonal_s$ari_onset - seasonal_s$v1_day
seasonal_s$v1_day <- seasonal_s$v1_day - symp_seasonal_s$ari_onset - seasonal_s$v1_day

# remove asymptomatics, subclinical and no onset
seasonal_s <- seasonal_s[seasonal_s$v1_day != 0, ]

plot_seasonal_s <- data.frame(hhID = rep(NA, (dim(seasonal_s)[1])*3))
plot_seasonal_s$hhID <- c(seasonal_s$hhID, seasonal_s$hhID, seasonal_s$hhID)
plot_seasonal_s$pcr <- c(seasonal_s$qPCR.1, seasonal_s$qPCR.2, seasonal_s$qPCR.3)
plot_seasonal_s$day <- c(seasonal_s$v1_day, seasonal_s$v2_day, seasonal_s$v3_day)
plot_seasonal_s$adult <- c(seasonal_s$adult, seasonal_s$adult, seasonal_s$adult)
plot_seasonal_s <- plot_seasonal_s[plot_seasonal_s$pcr > 0, ]
plot_seasonal_s$pcr <- log10(plot_seasonal_s$pcr)
plot_seasonal_s <- plot_seasonal_s[!is.na(plot_seasonal_s$hhID),]

# mean shedding by day
mean_seasonal_s <- data.frame(day = rep(NA, 14))
for(i in 1:14){
  mean_seasonal_s$day[i] <- (i-4)
  geo <- plot_seasonal_s$pcr[plot_seasonal_s$day==(i-4)]
  mean_seasonal_s$gmean[i] <- mean(geo)
}

# seasonal secondary symptoms
se_symp_s <- seasonal_s[c(1,2)]
se_symp_s <- merge(se_symp_s, symp)
se_symp_s <- se_symp_s[order(se_symp_s$hhID, se_symp_s$member, se_symp_s$day), ]

for (i in 1:dim(se_symp_s)[1]){
    se_symp_s$up[i] <- sum(se_symp_s$sthroat[i], se_symp_s$rnose[i], na.rm=TRUE)
    se_symp_s$lo[i] <- sum(se_symp_s$cough[i], se_symp_s$phlegm[i], na.rm=TRUE)
    se_symp_s$sy[i] <- sum(se_symp_s$fever[i], se_symp_s$headache[i], se_symp_s$pmuscle[i], na.rm=TRUE)
}
    
# add date of ARI onset
se_symp_s <- merge(se_symp_s, symp_seasonal_s[c(1,2,13)], by=c("hhID","member"))
se_symp_s$day <- se_symp_s$day - se_symp_s$ari_onset

# seasonal secondary symptom plot
se_s_symp_plot <- data.frame(day=rep(NA, 14))

for (i in 1:14) {
  se_s_symp_plot$day[i] <- (i-4)
  se_s_symp_plot$up[i] <- mean(se_symp_s$up[se_symp_s$day == (i-4)], na.rm=TRUE)
  se_s_symp_plot$lo[i] <- mean(se_symp_s$lo[se_symp_s$day == (i-4)], na.rm=TRUE)
  se_s_symp_plot$sy[i] <- mean(se_symp_s$sy[se_symp_s$day == (i-4)], na.rm=TRUE)
  se_s_symp_plot$temp[i] <- mean(se_symp_s$bodytemp[se_symp_s$day == (i-4)], na.rm=TRUE)
}

se_s_symp_plot$up <- se_s_symp_plot$up/2
se_s_symp_plot$lo <- se_s_symp_plot$lo/2
se_s_symp_plot$sy <- se_s_symp_plot$sy/3

### TCID seasonal influenza
# swine index tcid50
seas_i_tcid <- seasonal_i[c("hhID", "member", "adult", "v1_day", "v2_day", "v3_day")]
seas_i_tcid <- merge(seas_i_tcid, tcid, by=c("hhID", "member"))

# seasonal index tcid50 plot
sei_tcid_plot <- data.frame(day=rep(NA,10), tcid=rep(NA,10))

for (i in 1:10){
  sei_tcid_plot[i,1] <- i
  sei_tcid_plot[i,2] <- mean(c(seas_i_tcid$q_culture.1[seas_i_tcid$v1_day == (i)], seas_i_tcid$q_culture.2[seas_i_tcid$v2_day == (i)], seas_i_tcid$q_culture.3[seas_i_tcid$v3_day == (i)]), na.rm=TRUE)
}

# seasonal secondary tcid50
seas_s_tcid <- seasonal_s[c("hhID", "member", "adult", "v1_day", "v2_day", "v3_day")]
seas_s_tcid <- merge(seas_s_tcid, tcid, by=c("hhID", "member"))

# swine sceondary tcid50 plot
ses_tcid_plot <- data.frame(day=rep(NA,11), tcid=rep(NA,11))

for (i in 1:11){
  ses_tcid_plot[i,1] <- (i-4)
  ses_tcid_plot[i,2] <- mean(c(seas_s_tcid$q_culture.1[seas_s_tcid$v1_day == (i-4)], seas_s_tcid$q_culture.2[seas_s_tcid$v2_day == (i-4)], seas_s_tcid$q_culture.3[seas_s_tcid$v3_day == (i-4)]), na.rm=TRUE)
}


###                                  ###
###             Plot      (3X2)      ###
###                                  ###

windows(height=8, width=12)
layout(matrix(1:6, ncol=3, byrow=FALSE))
par(mar=c(1.5,3,1,1), oma=c(2,0,0,0))

# INDEX

# plot swine flu index
plot(0,0, axes=FALSE, xlim=c(0, 10), ylim=c(2, 8))
    lines(c(0, 10), c(log10(450), log10(450)), col="grey")
    lines(mean_i$day[c(2:10)], mean_i$gmean[c(2:10)], lwd=2)
    axis(1, pos=2, lwd=1.5, font=1, at=c(0, 2, 4, 6, 8, 10), label=c(0, 2, 4, 6, 8, 10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1, at=c(2:8),
         label=c(expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),
         expression(10^7),expression(10^8)),cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(plot_i$day[plot_i$adult == 1]), plot_i$pcr[plot_i$adult == 1], pch=1, cex=0.75)
    points(jitter(plot_i$day[plot_i$adult == 0]), plot_i$pcr[plot_i$adult == 0], pch=3, cex=0.75, col="red")
    mtext("Viral load by RT-PCR (copies/ml)", side = 2, line = 1.9, cex=0.7, font=2)  
    text(2.5,8.2,"Viral shedding by RT-PCR
                (pandemic)",adj=c(0,1),font=2, cex=1.1)                
    legend (8,8, c("child","adult"), pch=c(3,1), col=c("red","black"), cex=1.1, bty="n")

# plot seasonal flu index
plot(0,0, axes=FALSE, xlim=c(0, 10), ylim=c(2, 8))
    lines(c(0, 10), c(log10(450), log10(450)), col="grey")
    lines(mean_seasonal_i$day[c(2:10)], mean_seasonal_i$gmean[c(2:10)], lwd=2)
    axis(1, pos=2, lwd=1.5, font=1, at=c(0, 2, 4, 6, 8, 10), label=c(0, 2, 4, 6, 8, 10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1, at=c(2:8),
         label=c(expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6),
         expression(10^7),expression(10^8)),cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(plot_seasonal_i$day[plot_seasonal_i$adult == 1]), plot_seasonal_i$pcr[plot_seasonal_i$adult == 1], pch=1, cex=0.75)
    points(jitter(plot_seasonal_i$day[plot_seasonal_i$adult == 0]), plot_seasonal_i$pcr[plot_seasonal_i$adult == 0], pch=3, cex=0.75, col="red")
    mtext("Viral load by RT-PCR (copies/ml)", side = 2, line = 1.9, cex=0.7, font=2)  
    text(2.5,8.2,"Viral shedding by RT-PCR
                (seasonal)",adj=c(0,1),font=2, cex=1.1)

# plot swine flu index
plot(-3,0, axes=FALSE, xlim=c(0, 10), ylim=c(0, 5))
    lines(c(0, 10), c(0.3, 0.3), col="grey")
    lines(si_tcid_plot$day, si_tcid_plot$tcid, lwd=2)
    axis(1, pos=0, lwd=1.5, font=1, at=c(0, 2, 4, 6, 8, 10), label=c(0, 2, 4, 6, 8, 10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1, cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(c(swine_i_tcid$v1_day[swine_i_tcid$adult == 1], swine_i_tcid$v2_day[swine_i_tcid$adult == 1], swine_i_tcid$v3_day[swine_i_tcid$adult == 1]))
           , c(swine_i_tcid$q_culture.1[swine_i_tcid$adult == 1], swine_i_tcid$q_culture.2[swine_i_tcid$adult == 1], swine_i_tcid$q_culture.3[swine_i_tcid$adult == 1])
           , pch=1, cex=0.75)
    points(jitter(c(swine_i_tcid$v1_day[swine_i_tcid$adult == 0], swine_i_tcid$v2_day[swine_i_tcid$adult == 0 ], swine_i_tcid$v3_day[swine_i_tcid$adult == 0]))
           , c(swine_i_tcid$q_culture.1[swine_i_tcid$adult == 0], swine_i_tcid$q_culture.2[swine_i_tcid$adult == 0], swine_i_tcid$q_culture.3[swine_i_tcid$adult == 0])
           , pch=3, cex=0.75, col="red")
    mtext(expression(bold(paste("Viral load by culture (",log[10]," ",TCID[50],")", sep=""))), side = 2, line = 1.9, cex=0.7, font=2)  
    text(2.5,5.2,"Viral shedding by culture
               (pandemic)",adj=c(0,1),font=2, cex=1.1)    
    legend (8,5, c("child","adult"), pch=c(3,1), col=c("red","black"), cex=1.1, bty="n")

# plot seasonal flu index
plot(-3,0, axes=FALSE, xlim=c(0, 10), ylim=c(0, 5))
    lines(c(0, 10), c(0.3, 0.3), col="grey")
    lines(sei_tcid_plot$day, sei_tcid_plot$tcid, lwd=2)
    axis(1, pos=0, lwd=1.5, font=1, at=c(0, 2, 4, 6, 8, 10), label=c(0, 2, 4, 6, 8, 10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1, cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(c(seas_i_tcid$v1_day[seas_i_tcid$adult == 1], seas_i_tcid$v2_day[seas_i_tcid$adult == 1], seas_i_tcid$v3_day[seas_i_tcid$adult == 1]))
           , c(seas_i_tcid$q_culture.1[seas_i_tcid$adult == 1], seas_i_tcid$q_culture.2[seas_i_tcid$adult == 1], seas_i_tcid$q_culture.3[seas_i_tcid$adult == 1])
           , pch=1, cex=0.75)
    points(jitter(c(seas_i_tcid$v1_day[seas_i_tcid$adult == 0], seas_i_tcid$v2_day[seas_i_tcid$adult == 0], seas_i_tcid$v3_day[seas_i_tcid$adult == 0]))
           , c(seas_i_tcid$q_culture.1[seas_i_tcid$adult == 0], seas_i_tcid$q_culture.2[seas_i_tcid$adult == 0], seas_i_tcid$q_culture.3[seas_i_tcid$adult == 0])
           , pch=3, cex=0.75, col="red")
    mtext(expression(bold(paste("Viral load by culture (",log[10]," ",TCID[50],")", sep=""))), side = 2, line = 1.9, cex=0.7, font=2) 
    text(2.5,5.2,"Viral shedding by culture
               (seasonal)",adj=c(0,1),font=2, cex=1.1)

# plot swine flu index symptoms
plot(0,-1, axes=FALSE, ylim=c(0, 1), xlim=c(0, 10), bty="n", xlab="", ylab="")
    axis(1, pos=0, lwd=1.5, font=1, at=c(0,2,4,6,8,10), label=c(0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1,las=1, cex.axis=1.3, mgp=c(1.5,0.7,0)) 
    lines(sw_i_symp_plot$day[c(1:11)], sw_i_symp_plot$up[c(1:11)], type="l", col="red", lty="dashed", lwd=2)
    lines(sw_i_symp_plot$day[c(1:11)], sw_i_symp_plot$lo[c(1:11)], type="l", col="blue", lty="dotdash", lwd=2)
    lines(sw_i_symp_plot$day[c(1:11)], sw_i_symp_plot$sy[c(1:11)], type="l", col="black", lwd=2)
    mtext("Symptom scores", side = 2, line = 1.9, cex=0.7, font=2) 
    text(2.5,1.04,"Symptom scores
     (pandemic)",adj=c(0,1),font=2, cex=1.1)    
    legend (6.4,1.0, c("lower respir.","upper respir.","systemic"), lwd=1.7, lty=c("dotdash","dashed","solid"), 
            col=c("blue","red","black"), cex=0.9, bty="n") 
    
# plot seasonal flu index symptoms
plot(0,-1, axes=FALSE, ylim=c(0, 1), xlim=c(0, 10), bty="n", xlab="", ylab="")
    axis(1, pos=0, lwd=1.5, font=1, at=c(0,2,4,6,8,10), label=c(0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5, 0.6, 0))
    axis(2, pos=0, lwd=1.5, font=1,las=1, cex.axis=1.3, mgp=c(1.5,0.7,0)) 
    lines(se_i_symp_plot$day[c(1:11)], se_i_symp_plot$up[c(1:11)], type="l", col="red", lty="dashed", lwd=2)
    lines(se_i_symp_plot$day[c(1:11)], se_i_symp_plot$lo[c(1:11)], type="l", col="blue", lty="dotdash", lwd=2)
    lines(se_i_symp_plot$day[c(1:11)], se_i_symp_plot$sy[c(1:11)], type="l", col="black", lwd=2)
    mtext("Symptom scores", side = 2, line = 1.9, cex=0.7, font=2) 
    text(2.5,1.04,"Symptom scores
     (seasonal)",adj=c(0,1),font=2, cex=1.1)
            
mtext("Time (days since illness onset)", side = 1, at=0.17, line = -0.4, cex=0.7, font=2, outer=TRUE)
mtext("Time (days since illness onset)", side = 1, at=0.50, line = -0.4, cex=0.7, font=2, outer=TRUE)
mtext("Time (days since illness onset)", side = 1, at=0.85, line = -0.4, cex=0.7, font=2, outer=TRUE)

# End of script
