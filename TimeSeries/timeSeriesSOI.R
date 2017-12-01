library(marima)
library(ggplot2)
library(dplyr)

data_library(readr)
inflows <- read.csv(file="~/Dropbox/WORKSPACE/SDDP/TimeSeries/inflows4.txt", header=TRUE, sep="\t")
series<-t(inflows)
# Define marima model
Model5 <- define.model(kvar=4,ar=c(1),ma=0)
# Estimate marima model
Marima5 <- marima(series,Model5$ar.pattern)  
# Calculate residuals by filtering
Resid <- arma.filter(series, Marima5$ar.estimates)
# Compare residuals
plot(c(1:100), Resid$residuals[3, 1:100],
     xlab='marima residuals', ylab='arma.filter residuals')

ts.plot(data.frame(series[2,2:100], Marima5$fitted[2,2:100]) , gpars = list(col = c("black", "red")))
Marima5$ar.estimates
Marima5$Constant



library(MASS)
library(transport)
m = 4
num_outcomes = 3
fit = list()
rhsnoice =matrix(0, m, num_outcomes)

for(i in 1:m){
  fit[[i]]<-MASS::fitdistr(Marima5$residuals[i,2:103], densfun = "normal")
  rhsnoice[i,1] <- Marima5$Constant[i] + fit[[i]]$estimate[1] - fit[[i]]$estimate[2] 
  rhsnoice[i,2] <- Marima5$Constant[i] + fit[[i]]$estimate[1] 
  rhsnoice[i,3] <- Marima5$Constant[i] + fit[[i]]$estimate[1] + fit[[i]]$estimate[2] 
}


mm  = -Marima5$ar.estimates[,,2]
write.csv(mm,file = "./AR1Matrix.csv")
v0 = series[,1]

v1 = (mm%*%v0) + matrix(rhsnoice[,3],4,1)
v1 = (mm%*%v0) + Marima5$Constant
v1 = (mm%*%v1) + matrix(rhsnoice[,2],4,1)

#TODO: mirar si puedo hace los wasserstein 
f1 <- MASS::fitdistr(Marima5$residuals[i,2:103], densfun = "normal")
set.seed(27)
x <- pp(matrix(runif(500),250,2))
y <- pp(matrix(runif(500),250,2))
wasserstein(x,y,p=1)
wasserstein(x,y,p=2)



arsol <- ar(inflows[,2:5], method = "burg")
x=ts(inflows) #this makes sure R knows that x is a time series
x = x[,2:5]
plot(x, type="b") #time series plot of x with points marked as “o”
install.packages("astsa")
library(astsa) # See note 1 below
lag1.plot(x,1) # Plots x versus lag 1 of x.
acf(x, xlim=c(1,500)) # Plots the ACF of x for lags 1 to 19
xlag1=lag(x,-1) # Creates a lag 1 of x variable. See note 2
y=cbind(x,xlag1) # See note 3 below
ar1fit=lm(y[,1:4]~y[,5:8])#Does regression, stores results object named ar1fit
summary(ar1fit) # This lists the regression results
plot(ar1fit$fit,ar1fit$residuals) #plot of residuals versus fits
acf(ar1fit$residuals, xlim=c(1,18)) # ACF of the residuals for lags 1 to 18