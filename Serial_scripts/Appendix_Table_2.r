#
# R syntax to reproduce information for Figure 1 from:
#
# Cowling BJ, Fang VJ, Riley S, Peiris JS, Leung GM. 
# An estimate of the serial interval of influenza using 
# laboratory-confirmed natural infections in households. 
# Epidemiology, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Feburary 2, 2009

sedcase <- read.csv("http://web.hku.hk/~bcowling/influenza/Serial_scripts/serial_appendix_table_1.csv", header=TRUE)
source("http://www.hku.hk/bcowling/influenza/Serial_scripts/functions.r")

appt2 <- matrix(rep(NA,9),ncol=3,byrow=TRUE)

# Fit the weibull/gamma/lognormal model (allow for truncted serial interval)

serial.weibull.s <- optim(c(0.5,1.5), weibull.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
appt2[1,] <- c(round(exp(serial.weibull.s$par[2])*gamma(1+exp(-serial.weibull.s$par[1])),1),     # the mean
               round(exp(2*serial.weibull.s$par[2])*(gamma(1+2*exp(-serial.weibull.s$par[1]))
                 -(gamma(1+exp(-serial.weibull.s$par[1])))^2),1),                              # the variance
               round(serial.weibull.s$value*2+2*2,1))                                              # AIC

serial.gamma.s <- optim(c(0.5,0.5), gamma.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
appt2[2,] <- c(round(exp(serial.gamma.s$par[1])/exp(serial.gamma.s$par[2]),1),     # the mean
               round(exp(serial.gamma.s$par[1])/exp(2*serial.gamma.s$par[2]),1),   # the variance
               round(serial.gamma.s$value*2+2*2,1))                                  # AIC

serial.lognorm.s <- optim(c(0.5,-1.5), lognorm.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
appt2[3,] <- c(round(exp(exp(serial.lognorm.s$par[1])+exp(2*serial.lognorm.s$par[2])/2),1),      # the mean
               round(exp(exp(2*serial.lognorm.s$par[2])+2*exp(serial.lognorm.s$par[1]))
                                          *(exp(exp(2*serial.lognorm.s$par[2]))-1),1),         # the variance
               round(serial.lognorm.s$value*2+2*2,1))                                              # AIC

rownames(appt2) <- c("Weibull","Gamma","Lognormal")
colnames(appt2) <- c("Mean","Variance","AIC")
appt2

# End of script
