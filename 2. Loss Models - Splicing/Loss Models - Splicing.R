# This file is part of the Actuarial Statistics Project.
#
# Actuarial Statistics project is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Actuarial Statistics project is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Acturial Statistics Project  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015 Jellen Vermeir 
# jellenvermeir@gmail.com

# Set the working directory
# setwd("..")

library(MASS)
secura <- read.table(file="data/secura1.txt",header=TRUE)
attach(secura)
dim(secura)
names(secura)
head(secura)

############################################################################
################## DESCRIPTIVE STATISTICS ##################################
############################################################################

layout(matrix(c(1,3,2,4),nrow=2,ncol=2))
# ECDF
plot(ecdf(Loss), do.points = FALSE, xlab = 'Claim size', 
     ylab = 'ECDF', main = "Secura Belgian Re - Empirical PDF",col="blue")
hist(Loss,col="green",main = "Secura Belgian Re - Histogram of Losses")
plot(Year, Loss, ylab = "Claim size", xlab = "Year", 
     main = "Secura Belgian Re - Losses per year")
boxplot(Loss, main="Secura Belgian Re - Boxplot")

summary(Loss)

# lower truncation point
trunclower = 1200000


############################################################################
################## SINGLE DATA GENERATING PROCESS ##########################
############################################################################

# 1. Exponential likelihood
Exp.lik <- function(data,p)
{
  Exp.lambda <- exp(p[1])
  Exp.Lik <- dexp(data,rate=Exp.lambda,log=TRUE)
  return(-sum(Exp.Lik))
}

fit.nlm <- nlm(Exp.lik, p = c(log(1/mean(Loss-trunclower))),
               data=Loss-trunclower,print.level=2)
Exp.lambda <- exp(fit.nlm$estimate[1])

layout(matrix(c(1,3,2,3),nrow=2,ncol=2))
truehist(Loss, nbins = 40, ylim = c(0, 8e-07),xlab="Claim size",
         main = "Density estimate - Exponential")
curve(dexp(x-trunclower,rate=Exp.lambda),
      from = trunclower, to = 8e+06, n = 10000, 
      add = TRUE, col = "blue", lwd = 2)
legend('topright', legend = c("Exponential Density", 
            "Observed Relative Frequency"), col = c("blue", "cyan"), 
       pch=c(NA,15), pt.cex=2, lty = c(1, NA), lwd = c(2, NA))

plot(ecdf(Loss), do.points = FALSE, xlab = 'Claim size', ylab = 'CDF', 
     main = 'ECDF versus fitted CDF - Exponential', 
     xlim = c(trunclower, max(Loss)), lwd = 2)
curve((x >= trunclower) * pexp(x-trunclower,rate=Exp.lambda), 
      n = 1000, col = "blue", lwd = 2, add = TRUE)
legend('right', c('ECDF', 'Fitted CDF'), col = c(1, 4), lwd = 2)


qExp <- function(pinput,data){
  # solve F_X(x) = p numerically by minimizing (F_X(x) - p)^2
  objective <- function(xi){
    x <- exp(xi)
    F_X=0
    if(x >= 0)
      F_X <- pexp(x,rate=Exp.lambda)
    return((F_X-pinput)^2)
  }    
  optim(par = log(qexp(pinput, 1/mean(data))), 
        objective, method = "L-BFGS-B")$par
}

n <- length(Loss)
p <- 1:n/(n+1)
q <- exp(sapply(p,qExp,data=Loss-trunclower) )+ trunclower
q <- qexp(p,rate=Exp.lambda)+trunclower
plot(q, sort(Loss), 
     xlab = "Fitted quantiles", ylab = "Empirical quantiles", pch=1,
     main = "QQ plot - Exponential Distribution")
abline(a = 0, b = 1, lwd = 2, col = "red")

AIC.EXP <- (+2)*fit.nlm$minimum + 2*length(fit.nlm$estimate)
AIC.EXP


#2. lognormal likelihood
LN.lik <- function(data,p)
{
  LN.mu <- p[1]
  LN.sigma <- exp(p[2])
  LN.Lik <- dlnorm(data,meanlog=LN.mu,sdlog=LN.sigma,log=TRUE)
  # LN.Lik <- -log(data) - log(LN.sigma) -
  #                   0.5*log(2*pi)-0.5*(log(data)-LN.mu)^2/LN.sigma^2
          
  return(-sum(LN.Lik))
}

fit.nlm <- nlm(LN.lik, p = c(mean(log(Loss-trunclower)), 
              log(sd(log(Loss-trunclower)))),data=Loss-trunclower,
                                                    print.level=2)
LN.mu <- fit.nlm$estimate[1]
LN.sigma <- exp(fit.nlm$estimate[2])

truehist(Loss, nbins = 40, ylim = c(0, 8e-07),
         main = "Density estimate - Lognormal",
                                  xlab = "Claim size")
curve(dlnorm(x-trunclower,meanlog=LN.mu,sdlog=LN.sigma),
                main = "Density estimate - Lognormal",
                from = trunclower, to = 8e+06, n = 10000, 
                        add = TRUE, col = "blue", lwd = 2)
legend('topright', legend = c("Lognormal Density", 
          "Observed Relative Frequency"), col = c("blue", "cyan"), 
            pch=c(NA,15), pt.cex=2, lty = c(1, NA), lwd = c(2, NA))

plot(ecdf(Loss), do.points = FALSE, xlab = 'Claim size', ylab = 'CDF', 
     main = "ECDF versus fitted CDF - Lognormal", 
                      xlim = c(trunclower, max(Loss)), lwd = 2)
curve((x >= trunclower) * plnorm(x-trunclower,meanlog=LN.mu,sdlog=LN.sigma), 
      n = 1000, col = "blue", lwd = 2, add = TRUE)
legend('right', c('ECDF', 'Fitted CDF'), col = c(1, 4), lwd = 2)


qLN1 <- function(pinput,data){
  # solve F_X(x) = p numerically by minimizing (F_X(x) - p)^2
  objective <- function(xi){
    x <- exp(xi)
    F_X=0
    if(x >= 0)
      F_X <- plnorm(x,meanlog=LN.mu,sdlog=LN.sigma)
    return((F_X-pinput)^2)
  }    
  optim(par = log(qlnorm(pinput, mean(log(data)), sd(log(data)))), 
                        objective, method = "L-BFGS-B")$par
}

n <- length(Loss)
p <- 1:n/(n+1)
q <- exp(sapply(p,qLN1,data=Loss-trunclower)) + trunclower
q <- qlnorm(p,meanlog=LN.mu,sdlog=LN.sigma)+trunclower
plot(q, sort(Loss), main = "QQ plot - Lognormal distribution", 
     xlab = "Fitted quantiles", ylab = "Empirical quantiles", pch=1)
abline(a = 0, b = 1, lwd = 2, col = "red")


AIC.LN <- (+2)*fit.nlm$minimum + 2*length(fit.nlm$estimate)
AIC.LN

################################################################
#####################    BODY-TAIL #############################
################################################################

# observations equal to X_{n-k, n}
Loss[Loss == threshold]
# threshold specified using EVT technique 
threshold = 2580026
# determine splicing weights (using bernouilli MLE)
n = length(Loss)
k = length(Loss[Loss > threshold])
# body component
pn <- (n-k)/n
# tail component
k/n

# 3. Body-Tail: Exponential + GPD
EXPGDP.lik = function(p){
  EXP.lambda <- exp(p[1])
  
  GDP.shape <- p[2]
  # GDP.sigma <- exp(p[3])
  
  GDP.sigma <- (1-pn)*pexp(threshold-trunclower,rate=EXP.lambda)/
    (pn*dexp(threshold-trunclower,rate=EXP.lambda))
  
  seqThreshold <- which(Loss <= threshold)
  gtThreshold <- which(Loss > threshold)
  
  L1 <- pn * dexp(Loss[seqThreshold]-trunclower,
                                rate=EXP.lambda) / 
              pexp(threshold-trunclower, rate=EXP.lambda)
  L2 <- (1-pn) * (1/GDP.sigma)*(1 + (GDP.shape*
                          (Loss[gtThreshold]-threshold)) /
                            GDP.sigma)^(-1-1/GDP.shape)
  
  # Log Likelihood
  logL <- sum(c(log(c(L1,L2))))
  return(-logL)
}

# fit.nlm <- nlm(EXPGDP.lik, p = c(log(1/mean(Loss-trunclower)),
#                                  1/4,log(748329)),print.level=2)
fit.nlm <- nlm(EXPGDP.lik, p = c(log(1/mean(Loss-trunclower)),
                                 1/4),print.level=2)
Exp.lambda = exp(fit.nlm$estimate[1])
GDP.shape <- fit.nlm$estimate[2]
# GDP.sigma <- exp(fit.nlm$estimate[3])

GDP.sigma <- (1-pn)*pexp(threshold-trunclower,rate=Exp.lambda)/
  (pn*dexp(threshold-trunclower,rate=Exp.lambda))

layout(matrix(c(1,3,2,3),nrow=2,ncol=2))
truehist(Loss, nbins = 40, ylim = c(0, 8e-07),
         xlab="Claim size", main="Splicing - Continuation")
curve(I(x <= threshold) * pn * dexp(x-trunclower,rate=Exp.lambda) / 
        pexp(threshold-trunclower, rate=Exp.lambda) +
        (1-I(x <= threshold)) * (1-pn) * (1/GDP.sigma)*
        (1 + (GDP.shape*(x-threshold))/GDP.sigma)^(-1-1/GDP.shape), 
      from = trunclower, to = 8e+06, n = 10000, 
      add = TRUE, col = "blue", lwd = 2)
legend('topright', legend = c("Spliced Density", 
                              "Observed Relative Frequency"), 
              col = c("blue", "cyan"), pch=c(NA,15), pt.cex=2, 
                                lty = c(1, NA), lwd = c(2, NA))

plot(ecdf(Loss), do.points = FALSE, 
            xlab = 'Claim size', ylab = 'CDF', 
              main = 'Splicing - Exponential and GPD', 
                xlim = c(trunclower, max(Loss)), lwd = 2)
# Fitted CDF
curve((x >= trunclower) * ((x <= threshold) * pn * 
                pexp(x-trunclower,rate=Exp.lambda) / 
            pexp(threshold-trunclower, rate=Exp.lambda) + 
          (x > threshold)*(pn + (1-pn)*(1-(1+ GDP.shape *
                (x-threshold)/GDP.sigma)^(-1/GDP.shape)))), 
                  n = 1000, col = "blue", lwd = 2, add = TRUE)
legend('right', c('ECDF', 'Fitted CDF'), col = c(1, 4), lwd = 2)

## QQ plot
# this funtion takes a probablity p as argument 
# and returns the quantile x for which F_X(x) = p
qsplicing <- function(p){
  # solve F_X(x) = p numerically by minimizing (F_X(x) - p)^2
  objective <- function(xi){
    x <- exp(xi)
    F_X=500
    if(x >= trunclower & x <= threshold)
      F_X <- pn * pexp(x-trunclower,rate=Exp.lambda) / 
      pexp(threshold-trunclower, rate=Exp.lambda)
    if(x > threshold)
      F_X <- (pn + (1-pn)*
        (1-(1+ GDP.shape*(x-threshold)/GDP.sigma)^(-1/GDP.shape)))
    return((F_X-p)^2)
    
  }    
  if(p < pn)
    init <- log(qexp(p,rate=Exp.lambda)+trunclower)
  else
    init <- log(qlnorm(p,mean(log(Loss)),sd(log(Loss))))
  optim(par = init, objective, method = "L-BFGS-B")$par
}

n <- length(Loss)
p <- 1:n/(n+1)
q <- exp(sapply(p, qsplicing))
plot(q, sort(Loss), 
     main = "QQ plot - Splicing (Exponential and GPD)", 
     xlab = "Fitted quantiles", ylab = "Empirical quantiles", 
                                                        pch=1)
abline(a = 0, b = 1, lwd = 2, col = "red")

AIC.EXPGPD <- (+2)*fit.nlm$minimum + 2*length(fit.nlm$estimate)
AIC.EXPGPD

# 4. Body-Tail: Lognormal + GPD
LNGDP.lik = function(p){
  LN.mu <- p[1]
  LN.sigma <- exp(p[2])
  
  GDP.shape <- p[3]
  GDP.sigma <- (1-pn)*
      plnorm(threshold-trunclower,meanlog=LN.mu,sdlog=LN.sigma)/
    (pn*dlnorm(threshold-trunclower,meanlog=LN.mu,sdlog=LN.sigma))
  
  seqThreshold <- which(Loss <= threshold)
  gtThreshold <- which(Loss > threshold)
  
  L1 <- pn * dlnorm(Loss[seqThreshold]-
                      trunclower,meanlog=LN.mu,sdlog=LN.sigma) / 
        plnorm(threshold-trunclower, meanlog=LN.mu,sdlog=LN.sigma)
  L2 <- (1-pn) * (1/GDP.sigma)*
                  (1 + (GDP.shape*(Loss[gtThreshold]-threshold)) /
                                        GDP.sigma)^(-1-1/GDP.shape)
  
  # Log Likelihood
  logL <- sum(c(log(c(L1,L2))))
  return(-logL)
}

fit.nlm <- nlm(LNGDP.lik, p = c(mean(log(Loss)), 
                    log(sd(log(Loss))),0.1),print.level=2)
LN.mu = fit.nlm$estimate[1]
LN.sigma = exp(fit.nlm$estimate[2])
GDP.shape <- fit.nlm$estimate[3]
GDP.sigma <- (1-pn)*
    plnorm(threshold-trunclower,meanlog=LN.mu,sdlog=LN.sigma)/
  (pn*dlnorm(threshold-trunclower,meanlog=LN.mu,sdlog=LN.sigma))

layout(matrix(c(1,3,2,3),nrow=2,ncol=2))
truehist(Loss, nbins = 40, ylim = c(0, 8e-07),
         xlab = "Claim size", main="Splicing - Lognormal and GPD")
curve(I(x <= threshold) * pn * 
          dlnorm(x-trunclower,meanlog=LN.mu,sdlog=LN.sigma) / 
    plnorm(threshold-trunclower, meanlog=LN.mu,sdlog=LN.sigma) +
                  (1-I(x <= threshold)) * (1-pn) * (1/GDP.sigma)*
        (1 + (GDP.shape*(x-threshold))/GDP.sigma)^(-1-1/GDP.shape), 
                          from = trunclower, to = 8e+06, n = 10000, 
                                 add = TRUE, col = "blue", lwd = 2)
legend('topright', legend = c("Spliced Density", 
    "Observed Relative Frequency"), col = c("blue", "cyan"), 
        pch=c(NA,15), pt.cex=2, lty = c(1, NA), lwd = c(2, NA))


plot(ecdf(Loss), do.points = FALSE, xlab = 'Claim size', 
              ylab = 'CDF', main = "Splicing - Lognormal & GPD", 
                        xlim = c(trunclower, max(Loss)), lwd = 2)
curve((x >= trunclower) * ((x <= threshold) * pn * 
            plnorm(x-trunclower,meanlog=LN.mu,sdlog=LN.sigma) / 
        plnorm(threshold-trunclower, meanlog=LN.mu,sdlog=LN.sigma) + 
                    (x > threshold) * (pn + (1-pn)*(1-(1+ GDP.shape*
                          (x-threshold)/GDP.sigma)^(-1/GDP.shape)))), 
                        n = 1000, col = "blue", lwd = 2, add = TRUE)
legend('right', c('ECDF', 'Fitted CDF'), col = c(1, 4), lwd = 2)

# this funtion takes a probablity p as argument 
# and returns the quantile x 
# for which F_X(x) = p
qsplicing <- function(p){
  # solve F_X(x) = p numerically by minimizing (F_X(x) - p)^2
  objective <- function(xi){
    x <- exp(xi)
    F_X=50000
    if(x >= trunclower & x <= threshold)
      F_X <- pn * plnorm(x-trunclower,meanlog=LN.mu,sdlog=LN.sigma) / 
          plnorm(threshold-trunclower, meanlog=LN.mu,sdlog=LN.sigma)
    if(x > threshold)
      F_X <- (pn + (1-pn)*
        (1-(1+ GDP.shape*(x-threshold)/GDP.sigma)^(-1/GDP.shape)))
    return((F_X-p)^2)
  }   
  if(p < pn)
    init <- log(qlnorm(p, LN.mu, LN.sigma)+trunclower)
  else
    init <- log(qlnorm(p, mean(log(Loss)), sd(log(Loss))))
  optim(par = init, objective, method="L-BFGS-B")$par
}
n <- length(Loss)
p <- 1:n/(n+1)
q <- exp(sapply(p, qsplicing))
plot(q, sort(Loss), main = "QQ plot - Splicing (Lognormal and GPD)", 
     xlab = "Fitted quantiles", ylab = "Empirical quantiles", pch=1)
abline(a = 0, b = 1, lwd = 2, col = "red")


AIC.LNGPD <- (+2)*fit.nlm$minimum + 2*length(fit.nlm$estimate)
AIC.LNGPD

AIC.EXP
AIC.LN
AIC.EXPGPD
AIC.LNGPD