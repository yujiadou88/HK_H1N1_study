#
# R syntax to reproduce information for Figure 2 from:
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

# Fit the weibull model (allow for truncted serial interval)

serial.weibull <- optim(c(1,1), weibull.loglik, time=sedcase$serial_interval, data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

bt <- 1000
set.seed(1234567)
wei.pboot <- pboot.weibull(bt, sedcase$serial_interval, sedcase)    
quantile(wei.pboot$mean, c(0.5, 0.025, 0.975))          # 95% CI



############### Plot the graph ######################

windows(width=6, height=4)
par(mar=c(4.2,4,1,1))

plot(0, xlim=c(0,10), ylim=c(0,0.40), type="n", xlab="Serial interval (days)", ylab="Density", bty="n", axes=FALSE)


for (i in 1:100){
  b <- data.frame(times=0:200/20, des.prob=dweibull(0:200/20, rep(exp(wei.pboot$par1[i]),201), rep(exp(wei.pboot$par2[i]),201) ))
   lines(b$times, b$des.prob, lty=1, col="grey", lwd=1)
}

b <- data.frame(times=0:200/20, des.prob=dweibull(0:200/20, rep(exp(serial.weibull$par[1]),201),rep(exp(serial.weibull$par[2]),201) ))
 lines(b$times, b$des.prob, lty=1, col="black", lwd=1.5)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)

# End of script
