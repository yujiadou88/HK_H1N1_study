#
# R syntax including functions needed in the analyses described in:
#
# Cowling BJ, Fang VJ, Riley S, Peiris JS, Leung GM. 
# An estimate of the serial interval of influenza using 
# laboratory-confirmed natural infections in households. 
# Epidemiology, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Feburary 2, 2009


# log likelihood functions ...

gamma.loglik <- function(parameters, time,data){               
    if(length(parameters)!=2) stop("gamma distribution should have two parameters")
    k <- exp(parameters[1])
    lambda <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dgamma(time.g[i], shape=k,scale=1/lambda) / (1-pgamma(trunct[i], shape=k,scale=1/lambda)) )
    }
    -sum(log.lik)
}


lognorm.loglik <- function(parameters, time, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    mean <- exp(parameters[1])
    sd <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dlnorm(time.g[i], meanlog=mean, sdlog=sd) / (1-plnorm(trunct[i],meanlog=mean, sdlog=sd)) )
    }
    -sum(log.lik)
}


weibull.loglik <- function(parameters, time, data){
    if(length(parameters)!=2) stop("lognormal distribution should have two parameters")
    alpha <- exp(parameters[1])
    beta <- exp(parameters[2])

    time.g <- na.exclude(time)
    n <- length(time.g)  
    trunct <- data$clinic_day[!is.na(time)]
    log.lik <- rep(NA,n)
    for(i in 1:n){
        log.lik[i] <- log( dweibull(time.g[i], shape=alpha, scale=beta) / (1-pweibull(trunct[i],shape=alpha, scale=beta)) )
    }
    -sum(log.lik)
}


### parametric bootstrap

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
    new.data <- new.data[new.data$day_symp>=new.data$clinic_day,]
    new.data <- new.data[1:n,]
    
    new.weibull <- optim(c(0.5,1.5), weibull.loglik, time=new.data$day_symp,data=new.data,
                    method="L-BFGS-B", lower=c(-Inf,-Inf), upper=c(Inf,Inf))
    output$mean[i] <- exp(new.weibull$par[2])*gamma(1+exp(-new.weibull$par[1]))
    output$par1[i] <- new.weibull$par[1]
    output$par2[i] <- new.weibull$par[2]
  }
  output
}



#
# Turnbull estimation involving interval censoring and left truncation
#

make.equivalence <- function(data)
{
    #
    # input data.frame should include columns "timeL" and "timeR"
    # n.b. in the code below, 0=left and 1=right
    # want to keep the last of a sequence of 0s, and the first of a sequence of 1s
    #
    m <- length(data[, 1])
    timeL <- matrix(c(data$timeL, rep(0, m)), byrow = F, ncol = 2)
    timeR <- matrix(c(data$timeR, rep(1, m)), byrow = F, ncol = 2)
    timeLR <- matrix(c(data$timeL, data$timeR, rep(0, m), rep(1, m), rep(0, 2 * m)),
        byrow = F, ncol = 3)
    newtimeLR <- timeLR[order(timeLR[, 1]),  ]
    for(i in 2:(2 * m)) {
        if(newtimeLR[i, 2] == 1 && newtimeLR[i - 1, 2]==1)
            newtimeLR[i, 3] <- 1
        if(newtimeLR[i, 2] == 0 && newtimeLR[i - 1, 2]==0)
            newtimeLR[i-1, 3] <- 1

    }
    temp.eq <- newtimeLR[(newtimeLR[, 3] != 1),  ]
    k <- length(temp.eq[, 1])/2
    eq.classes <- data.frame(timeL = temp.eq[2 * (1:k) - 1, 1], timeR = temp.eq[2 *
        (1:k), 1])
    eq.classes
}

make.eq.matrix <- function(data, eq.classes)
{
    #
    # takes data.frame "data" (containing columns timeL and timeR)
    # and mx2 matrix of associated equivalence classes "eq.classes"
    #
    n <- length(data[, 1])
    m <- length(eq.classes[, 1])
    output <- matrix(rep(NA, 2 * n * m ), nrow = 2*n)
    for(i in 1:n) {
        for(j in 1:m) {
            if((data$timeL[i] <= eq.classes[j, 1]) && (data$timeR[i] >= eq.classes[j, 2])) {
                output[i, j] <- 1
            }
            else {
                output[i, j] <- 0
            }
            if((data$trunct[i] <= eq.classes[j, 1]) ) {       # left trunctated, like index clinic visit day
                output[n+i, j] <- 1
            }
            else {
                output[n+i, j] <- 0
            }
        }
    }
    output
}

icnpF <- function(data, eq.class, eq.class.mat, EPS = 0.0001)
{
    # function to find non-parametric estimate of F (Turnbull, 1973)
    # takes data.frame "data" (containing columns timeL and timeR)
    # and mx2 matrix of associated equivalence classes "eq.classes"
    # and nxm matrix of indicators
    #
    n <- length(data[, 1])
    m <- length(eq.class[, 1])
    A <- eq.class.mat           # alpha_ij, beta_ij
    if(dim(A)[1] != n*2)
        stop("indicator matrix should have ", n*2, " rows!")
    if(dim(A)[2] != m)
        stop("indicator matrix should have ", m, " columns!")
    P <- matrix(rep(1/m, m), ncol = 1)          # s_j - set initial estimates s_j=1/m for all j
    iter <- 1
    Q <- matrix(rep(1, m), ncol = 1)
    mu <- matrix(rep(NA, n*m),ncol=m)
    v <- matrix(rep(NA, n*m),ncol=m)
    while(max(abs(Q - P)) > EPS) {
        iter <- iter + 1
        Q <- P
    for (j in 1:m){
         for (i in 1:n){
                 mu[i,j] <- A[i,j]*P[j] / A[i,] %*% P
         v[i,j] <- (1-A[n+i,j])*P[j] / A[n+i,] %*% P
         } 
        }
    M <- sum(mu+v)
    for (j in 1:m){
            P[j] <- sum(mu[,j]+v[,j])/M
    }
    }
    cat(iter, "iterations", "\n")
    data.frame(pi=P, cumpi=cumsum(P))
}


#
#
#

turnbull.lik <- function(pars, eq.mat){
    #
    # pars should have same dimension 
    #
    n <- dim(eq.mat)[1]/2
    m <- dim(eq.mat)[2]
    if(length(pars)!=m) stop("not the correct number of parameters")
    par.mat <- matrix(pars, ncol=1)
    log.lik <- log(eq.mat[1:n,] %*% par.mat) - log(eq.mat[1:n+n,] %*% par.mat) # + 10000*(sum(pars)-1)
    -1*sum(log.lik)
}

turnbull.lik.d1 <- function(pars, eq.mat){
    #
    # pars should have same dimension 
    #
    n <- dim(eq.mat)[1]/2
    m <- dim(eq.mat)[2]
    if(length(pars)!=m) stop("not the correct number of parameters")
    par.mat <- matrix(pars, ncol=1)
    log.lik.d1 <- rep(NA, m)
    for(j in 1:m){
        log.lik.d1[j] <- sum(eq.mat[1:n,j] / (eq.mat[1:n,] %*% par.mat)) - sum(eq.mat[1:n+n,j] / (eq.mat[1:n+n,] %*% par.mat))
    }
    log.lik.d1
}

turnbull.lik.d2 <- function(pars, eq.mat){
    #
    # pars should have same dimension 
    #
    n <- dim(eq.mat)[1]/2
    m <- dim(eq.mat)[2]
    if(length(pars)!=m) stop("not the correct number of parameters")
    par.mat <- matrix(pars, ncol=1)
    log.lik.d2 <- matrix(rep(NA, m*m), ncol=m)
    for(j in 1:m){
        for(k in 1:m){
            log.lik.d2[j,k] <- sum(-1*eq.mat[1:n,j]*eq.mat[1:n,k] / ( (eq.mat[1:n,] %*% par.mat) * (eq.mat[1:n,] %*% par.mat) ) ) 
                           - sum(-1*eq.mat[1:n+n,j]*eq.mat[1:n+n,k] / ( (eq.mat[1:n+n,] %*% par.mat) * (eq.mat[1:n+n,] %*% par.mat) ) ) 
        }
    }
    log.lik.d2
}

#
#
#

make.turnbull <- function(data){
    eqclass <- make.equivalence(data)
    m <- length(eqclass[,1])
    eqmat <- make.eq.matrix(data, eqclass)
    icnp.F <- icnpF(data, eqclass, eqmat)
    varcov <- -1*solve(turnbull.lik.d2(icnp.F[,1], eqmat))
    pi.se <- sqrt(diag(varcov))
    pi.lower <- rep(NA, m)
    pi.upper <- rep(NA, m)
    for(j in 1:m){
        pi.lower[j] <- max(0, icnp.F[j,1]-1.96*pi.se[j])
        pi.upper[j] <- min(1, icnp.F[j,1]+1.96*pi.se[j])
    }

    cum.pi.se <- rep(NA, m)
    for(j in 1:m){
        cum.pi.se[j] <- sqrt(sum(varcov[1:j,1:j]))
    }
    cum.pi.lower <- rep(NA, m)
    cum.pi.upper <- rep(NA, m)
    for(j in 1:m){
        cum.pi.lower[j] <- max(0, icnp.F[j,2]-1.96*cum.pi.se[j])
        cum.pi.upper[j] <- min(1, icnp.F[j,2]+1.96*cum.pi.se[j])
    }
    data.frame(eq.timeL=eqclass[,1], eq.timeR=eqclass[,2], pi=round(icnp.F[,1],6), 
        pi.se=round(pi.se,6), pi.lower=round(pi.lower,6), pi.upper=round(pi.upper,6),
    cum.pi=round(icnp.F[,2],6), cum.pi.se=round(cum.pi.se,6),
        cum.pi.lower=round(cum.pi.lower,6), cum.pi.upper=round(cum.pi.upper,6))
}

#
# plot survival function
#

turnbull.plot <- function(turnbull.object, line.type=1, line.col=1, line.width=1, max.time=100, CI.col=gray(0.8))
{
    #
    # function to add horizontal lines
    # vertical lines will be added only for equivalence classes with length 0
    #
    obj <- turnbull.object
    m <- length(obj$cum.pi)
        lines(x=c(0, obj$eq.timeL[1]), y=rep(0, 2), lty=line.type, col=line.col, lwd=line.width)
        for(i in 1:(m-1)){
            lines(c(obj$eq.timeR[i], obj$eq.timeL[i+1]), rep(obj$cum.pi[i], 2), lty=line.type,
                col=line.col, lwd=line.width)
        }
        lines(x=c(obj$eq.timeL[m], max.time), y=rep(1, 2), lty=line.type, col=line.col, lwd=line.width)
        if(obj$eq.timeL[1]==obj$eq.timeR[1]) lines(rep(obj$eq.timeL[1], 2),
            c(0, obj$cum.pi[1]), lty=line.type, col=line.col, lwd=line.width)
        for(i in 2:m){
            if(obj$eq.timeL[i]==obj$eq.timeR[i]) lines(rep(obj$eq.timeL[i], 2),
                c(obj$cum.pi[i-1], obj$cum.pi[i]), lty=line.type,
                    col=line.col, lwd=line.width)
        }
    invisible(0)
}


# End of script
