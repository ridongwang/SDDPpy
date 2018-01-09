library(marima)
library(ggplot2)
library(dplyr)
library(MASS)
library(transport)
setwd("~/Dropbox/WORKSPACE/SDDP/TimeSeries")

data_library(readr)

m = 50
inflows <- read.csv(file="~/Dropbox/WORKSPACE/SDDP/TimeSeries/inflows50.txt", header=TRUE, sep="\t")
series<-t(inflows)
# Define marima model
Model5 <- define.model(kvar=m,ar=c(1),ma=0)
# Estimate marima model
Marima5 <- marima(series,Model5$ar.pattern)  
# Calculate residuals by filtering
Resid <- arma.filter(series, Marima5$ar.estimates)
# Compare residuals
plot(c(1:100), Resid$residuals[3, 1:100],
     xlab='marima residuals', ylab='arma.filter residuals')

ts.plot(data.frame(series[30,2:1000], Marima5$fitted[30,2:1000]) , gpars = list(col = c("black", "red")))
Marima5$ar.estimates
Marima5$Constant


num_outcomes = 30
fit = list()
rhsnoice =matrix(0, m, num_outcomes)
step_size = 6.0/(num_outcomes-1)
sds = seq(from=-3, to=3, by=step_size)
for(i in 1:m){
  fit[[i]]<-MASS::fitdistr(Marima5$residuals[i,2:1000], densfun = "normal")
  for( j in 1:num_outcomes){
    rhsnoice[i,j] <- Marima5$Constant[i] + fit[[i]]$estimate[1] + sds[j]*fit[[i]]$estimate[2]
  }
}
write.csv(rhsnoice,file = "./RHSnoise50.csv")

mm  = -Marima5$ar.estimates[,,2]
write.csv(mm,file = "./AR1Matrix50.csv")



v0 = series[,1]

v1 = (mm%*%v0) + matrix(rhsnoice[,3],4,1)
v1 = (mm%*%v0) + Marima5$Constant
v1 = (mm%*%v1) + matrix(rhsnoice[,2],4,1)

####################
# Processing LP data
###################
library(ggplot2)
library(plotly)
setwd("~/Projects/SDDPOutputFiles")

temp = list.files(pattern="*.csv")
mdf = data.frame(matrix(ncol = 4, nrow = 0))
names(mdf)<- c("pass", "num_ctr", "lp_time", 'Algo')
for (i in 1:length(temp)) mdf<-rbind(mdf, read.csv(temp[i]))


model <- lm(mdf$lp_time~ poly(mdf$num_ctr,1))
summary(model)

plot(fitted(model),residuals(model))




gg=ggplot(model,aes(x=,y=price,color=color))
ggplotly(gg)

