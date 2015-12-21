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


xAxis = "Inter arrival time (days)"; mainTitle = "Inter arrival times"# Q1.
startDate <- as.Date("01/01/1980",format="%d/%m/%Y")
arrivalDates <- as.Date(read.table(file="data/DanishData.txt",
                                   header=TRUE)$Date, format="%m/%d/%Y")

# Arrival times since start of the observation period, expressed in days
Ti <- as.numeric(difftime(arrivalDates,startDate,units="days"))
# Interarrival times, expressed in days
Wi <- diff(c(0,Ti))
# Number of claims at time T (1 claim for each arrival time)
claimCount <- rep(1,length(Ti))

data <- data.frame(time=arrivalDates,Ti=Ti,Wi=Wi,
                    count=claimCount,year=format(arrivalDates,"%Y"))

#layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,1))
xAxis = "Time (days)"; yAxis = "cumulative # claims" 
mainTitle = "Cumulative number of arrivals"
n <- 1:length(claimCount)
plot(Ti,n,xlab=xAxis,ylab=yAxis,
     main=mainTitle)

hist(data$Wi,breaks=0:25-1,xlab=xAxis,main=paste(mainTitle,"(overall)"))

# Proposition.. homogeneous
# Add claimcounts that fall on same day together
# We now have rows that represent time interval N(t+s)-N(s) --> Wi
# and corresponding number of claims n --> count
library(plyr)
claims <- ddply(data[,c("time","Wi","count")],"time",numcolwise(sum))
claims$year <- format(claims$time,"%Y")

# Use the data to fit a poisson process
# Maximum likelihood for P[N(t+s)-N(s)=n]
Poi.lik = function(x,beta,expo){
  lambda = exp(beta) # transform unconstrained. Make lambda positive
  dens = exp(-lambda*expo)*(lambda*expo)^x/gamma(x+1)
  logl = sum(log(dens))
  return(-logl)
}
Poi.result = nlm(Poi.lik,1,hessian=TRUE,x=claims$count,expo=claims$Wi)
# Alternative
fm_pois <- glm(claims$count~1,family=poisson(link="log"),offset=log(claims$Wi)) 

lambda <- exp(Poi.result$estimate)
lambda
Poi.AIC <- (+2)*Poi.result$minimum+2*length(Poi.result$estimate) # -2*ln(L) + 2*k
Poi.AIC

lambda <- exp(fm_pois$coefficients[[1]])
lambda

# Or..
lambda = weighted.mean(claims$count/claims$Wi,claims$Wi)
lambda

# Or..
lambda = 1/mean(data$Wi)

par(mfrow=c(1,2))
simulation <- rexp(10000,rate=lambda)
endBreaks <- max(c(simulation,data$Wi))
h <- hist(data$Wi,breaks=0:ceiling(endBreaks+1)-0.5,plot=FALSE)
h$counts <- h$counts/nrow(data)
hSim <- hist(simulation,breaks=0:ceiling(endBreaks+1)-0.5,plot=FALSE)
hSim$counts <- hSim$counts/length(simulation)
plot(h, col=rgb(1,0,0,0.5),
     main="Frequency distribution of inter-arrivals",
     xlab="Claimcount", ylab="Normalized frequency")
lines(hSim, col=rgb(0,0,1,0.5))
legend("topright",lty=c(1,1),legend=c("Emp","Exp(0.54)"),col=c("blue","red"))

#boxplot(data$Wi, main="Boxplot of inter-arrivals")
mTitle="CDF of Inter-arrivals"
xAxis="Days between events";leg="Emp"
plot(ecdf(data$Wi),xlab=xAxis,ylab="CDF",
     col="blue",main=mTitle,xlim=c(-1,10))
lines(ecdf(rexp(10000,rate=lambda)),col="red")
legend("topright",lty=c(1,1),col=c("blue","red"),
       legend=c(leg,paste("Exp(",round(lambda,2),")",sep="")))

# Visualisation of inter arrival times, per year
par(mfrow=c(4,3))
xAxis = "Inter arrival time (days)"; mainTitle = "Inter arrival times"
hist(data$Wi,breaks=0:25-1,xlab=xAxis,main=paste(mainTitle,"(overall)"))
for(i in unique(data$year))
{
  indices <- which(data$year==i)
  hist(data$Wi[indices],breaks=0:25-1,
       xlab=xAxis,main=paste(mainTitle," (",i,")",sep=""))
}

# Cumulative claims, per year
par(mfrow=c(4,3))
xAxis = "Time (days)"; yAxis = "cumulative # claims" 
mainTitle = "Cumulative Arrivals"
n <- 1:length(claimCount)
plot(Ti,n,xlab=paste(xAxis,"(all years)"),ylab=yAxis,
     main=paste(mainTitle,"(All years)"))
for(i in unique(data$year))
{
  indices <- which(data$year==i)
  plot(Ti[indices],n[indices],xlab=paste(xAxis," (",i,")",sep=""),
       ylab=yAxis,main=paste(mainTitle," (",i,")",sep=""))
}

summaryStats <- data.frame(year=character(),
                 NrSamples=double(),
                 Minimum=double(),
                 First_Quartile=double(),
                 Median=double(),
                 Mean=double(),
                 Third_Quartile=double(),
                 Maximum=double(),
                 lambda=double(),
                 stringsAsFactors=FALSE)

summaryStats[1,] <- c("Overall",nrow(data),
                        summary(data$Wi),round(lambda,6))

par(mfrow=c(4,3))
plot(ecdf(data$Wi),xlab=xAxis,ylab="CDF",
     col="red",main=paste(mTitle,"(Overall)"),xlim=c(-1,10))
lines(ecdf(rexp(10000,rate=lambda)),col="blue")
legend("topright",lty=c(1,1),col=c("red","blue"),
       legend=c(leg,paste("Exp(",round(lambda,2),")",sep="")))
years <- unique(data$year)
for(i in 1:length(years))
{
  theYear <- years[[i]]
  indices <- which(data$year==theYear)
  lambda <- 1/mean(data$Wi[indices])
  
  summaryStats[i+1,] <- c(as.character(theYear),length(indices),
                          summary(data$Wi[indices]),round(lambda,6))
  
  plot(ecdf(data$Wi[indices]),xlab=xAxis,ylab="CDF",xlim=c(-1,10),
       col="red",main=paste(mTitle," (",theYear,")",sep=""))
  lines(ecdf(rexp(10000,rate=lambda)),col="blue")
  legend("topright",lty=c(1,1),col=c("red","blue"),
         legend=c(leg,paste("Exp(",round(lambda,2),")",sep="")))
}

AnnualMeans <- summaryStats$Mean
summaryStats

# install.packages("gridExtra")
library(gridExtra)
png("summaryStats.png", height=300, width=700)
p<-tableGrob(summaryStats)
grid.arrange(p)
dev.off()

# Create boxplot
par(mfrow=c(1,2))
boxplot(data$Wi ~ data$year,main="Distribution Inter-arrival times",
        xlab="Year",ylab="Inter-arrival times (days)")
# Number of claims per year
claimsPerYear <- ddply(data[,c("year","count")],"year",numcolwise(sum))
data <- data.frame(time=arrivalDates,Ti=Ti,Wi=Wi,count=claimCount)
data$year <- format(data$time,"%Y")
barplot(claimsPerYear$count,names.arg=claimsPerYear$year,
        main="Number of claims",xlab="Year",ylab="# Claims")

centeredMean <- function(center=51,width=50,data)
{
  if(!is.numeric(data))
    stop("Input data should be numeric vector")
  
  sum <- 0;
  n <- length(data)
  
  minIndex <- center-width
  maxIndex <- center+width
  
  sum = sum(data[max(1,minIndex):min(n,maxIndex)])
  
  # Account for proportion of indices below 1
  if(minIndex <= 0)
    sum = sum + (abs(minIndex)+1)*data[1]
  # Account for proportion of indices above n
  if(maxIndex > n)
    sum = sum + (maxIndex-n)*data[n]
  
  return(sum/(2*width+1))
}

m <- 50
Wi <- data$Wi
movingMeanInterArrival <- sapply(1:length(Wi),
                                 centeredMean,width=m,data=Wi)
movingLambda <- 1/movingMeanInterArrival

par(mfrow=c(2,1))
plot(data$time,movingMeanInterArrival,type="l",xlab="Time",
                        ylab="Mean Inter-Arrival",col="blue",
     main="Estimate of mean inter-arrival times")
plot(data$time,movingLambda,type="l",xlab="Time",
          ylab="Intensity",col="red",main="Estimate of intensity")

# Q5.
allDays <- seq(from=as.Date("1980/01/01"),to=as.Date("1990/12/31"),by="days")
lambdas <- rep(0,length(allDays));
for(i in unique(data$year))
{
  lambda <- as.numeric(summaryStats$lambda[summaryStats$year==i])
  indices <- which(format(allDays,"%Y")==i)
  lambdas[indices] <- lambda
}

t <- 1:length(lambdas)
mu_t <- sapply(t, function(x) sum(lambdas[1:x]))

par(mfrow=c(1,1))
plot(t,mu_t,xlab="Time t",type="l",col="blue",
     main=expression(paste("Continuous mean value function ",mu,"(0,t)",sep="")))


# Q6 and Q7
mu_Ti <- mu_t[Ti]
mu_Wi <- diff(c(0,mu_Ti))
mu_MovingMeanInterArrival <- sapply(1:length(mu_Wi),centeredMean,width=m,data=mu_Wi)
mu_MovingLambda <- 1/mu_MovingMeanInterArrival

empMean <- mean(mu_Wi)
empLambda <- 1/empMean
empLambda

# Lambda seems stationary around 1
par(mfrow=c(2,1))
plot(data$time,mu_MovingMeanInterArrival,type="l",xlab="Time",
     ylab="Mean Inter-Arrival",col="blue",
     main="Estimate of transformed mean inter-arrival times")
plot(data$time,mu_MovingLambda,type="l",xlab="Time",ylab="Intensity",
     col="red",main="Estimate of transformed intensity")

par(mfrow=c(2,2))
xAxis="Days between events";leg="Emp"
plot(ecdf(mu_Wi),xlab=xAxis,ylab="CDF",
     col="red",main=expression(paste("CDF of inter-arrivals"," ",
                                     mu,"(Ti)",sep="")),xlim=c(-1,10))
lines(ecdf(rexp(10000,rate=1)),col="blue")
legend("topright",lty=c(1,1),col=c("red","blue"),
       legend=c(leg,paste("Exp(1)",sep="")))

p <- ppoints(100)    # 100 equally spaced points on (0,1), excluding endpoints
q <- quantile(mu_Wi,p=p) # percentiles of the sample distribution
plot(qexp(p) ,q, main="Exponential Q-Q Plot",
     xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqline(q, distribution=qexp,col="blue", lty=2)

acf(mu_Wi, main = "ACF of inter-arrival times", lag.max = 20,
    ylab = "", xlab = "", col = "blue", ci.col = "red")
pacf(mu_Wi, main = "PACF of inter-arrival times", lag.max = 20,
     ylab = "", xlab = "", col = "blue", ci.col = "red")


lambda1980 <- as.numeric(summaryStats$lambda[summaryStats$year==1980])
data1980 <- data[data$year==1980,]

par(mfrow=c(2,2))
sumClaims <- cumsum(data1980$Wi)
Nt <- 0:length(sumClaims)
plot(stepfun(sumClaims,Nt),col="blue",xlab="Time (days)",
     ylab="Number of Claims",main="Paths of a Poisson Process (1980)")


# Simulate n interarrival times
paths <- 6; n <- nrow(data1980)
simulations <- matrix(rexp(paths*n, rate = lambda1980), nrow = paths)
t <- t(apply(simulations, 1, cumsum))
for(i in 1:paths)
{
  lines(stepfun(t[i,],Nt), xlab = "Time", 
                  ylab = "Number of claims", 
        main = "Paths of a Poisson Process",col="red")
}
legend("topleft",lty=c(1,1),col=c("blue","red"),
       legend=c("N(t) - 1980","Simulations"))


simulation <- rexp(10000,rate=lambda1980)

#par(mfrow=c(3,1))
endBreaks <- max(c(simulation,data1980$Wi))
h <- hist(data1980$Wi,breaks=0:ceiling(endBreaks+1)-0.5,plot=FALSE)
h$counts <- h$counts/nrow(data1980)
hSim <- hist(simulation,breaks=0:ceiling(endBreaks+1)-0.5,plot=FALSE)
hSim$counts <- hSim$counts/length(simulation)
plot(h, col=rgb(1,0,0,0.5),
     main="Inter-arrival times 1980 VS simulation",
     xlab="Inter-arrival times (Days)", ylab="Normalized frequency")
lines(hSim, col=rgb(0,0,1,0.5))
legend("topright",lty=c(1,1),
       legend=c("Inter-arrival times","Exp(0.455)"),col=c("blue","red"))

plot(density(data1980$Wi),xlim=c(-1,15),xlab="Inter-arrival times (Days)",
      ylim=c(0,0.35),col="blue",main="Inter-arrival times 1980 VS simulation")
lines(density(simulation),
     xlim=c(-1,15),ylim=c(0,0.35),col="red")
legend("topright",lty=c(1,1),
       legend=c("Inter-arrival times","Exp(0.455)"),col=c("blue","red"))

p <- ppoints(100)    # 100 equally spaced points on (0,1), excluding endpoints
q <- quantile(data1980$Wi,p=p) # percentiles of the sample distribution
plot(qexp(p,rate=lambda1980) ,q, main="Exponential Q-Q Plot",
     xlab=expression(paste("Theoretical Quantiles (",lambda,"=",0.455,")",sep="")),
                                                            ylab="Sample Quantiles")
qqline(q, distribution = function(p) qexp(p,rate=lambda1980),col="blue", lty=2)