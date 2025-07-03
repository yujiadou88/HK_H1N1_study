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

data.s <- sedcase[,1:2]
names(data.s)[2] <- "trunct"
data.s$timeL <- sedcase$serial_interval
data.s$timeR <- sedcase$serial_interval
data.s$exact <- 1
data.s <- data.s[!is.na(data.s$timeL),]
s.turnbull <- make.turnbull(data.s)

# Fit the weibull/gamma/lognormal model (allow for truncted serial interval)

serial.weibull.s <- optim(c(0.5,1.5), weibull.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

serial.gamma.s <- optim(c(0.5,0.5), gamma.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

serial.lognorm.s <- optim(c(0.5,-1.5), lognorm.loglik, time=sedcase$serial_interval,data=sedcase,
 method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))

# Plot the figure

windows(width=6, height=4)
par(mar=c(4,4,1,1))

plot(0, xlim=c(0,10), ylim=c(0,1), type="n", xlab="", ylab="", bty="n", axes=FALSE)

turnbull.plot(s.turnbull, line.type=3, line.col="black", line.width=1, max.time=10)

b <- data.frame(times=0:100/10, 
                cum.prob=plnorm(0:100/10, rep(exp(serial.lognorm.s$par[1]),101), rep(exp(serial.lognorm.s$par[2]),101)))
lines(b$times, b$cum.prob, lty=4, col="black", lwd=1)

b <- data.frame(times=0:100/10, 
                cum.prob=pgamma(0:100/10,rep(exp(serial.gamma.s$par[1]),101), rep(exp(serial.gamma.s$par[2]),101) ))
lines(b$times, b$cum.prob, lty=2, col="black", lwd=1)

b <- data.frame(times=0:100/10, 
                cum.prob=pweibull(0:100/10, rep(exp(serial.weibull.s$par[1]),101), rep(exp(serial.weibull.s$par[2]),101) ))
lines(b$times, b$cum.prob, lty=1, col="black", lwd=1)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=2*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)

mtext("Serial interval (days)", side=1, line=2.5, cex=1.0)  
mtext("Cumulative proportion", side=2, line=3, cex=1.0) 

# End of script
