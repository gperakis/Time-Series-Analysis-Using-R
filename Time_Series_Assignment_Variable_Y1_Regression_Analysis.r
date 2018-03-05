
# loads the library urca, which contains the command ur.df() 
# and applies Unit-Root Testing based on Dickey-Fuller test 
# install.packages(c('fGarch', 'urca','tseries','MASS'))
library(urca)
library(tseries)
library(MASS)
library(fGarch)

# Importing data into R can be carried out in various ways. Below, the command read.table is used:
data_all = read.table('/home/george/AUEB/TimeSeries/time_series_in_r/Data-Assignment-Feb-2018.txt')

xdata_cols <- c('V5','V6', 'V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18')

xdata <- data_all[c(xdata_cols)]

head(xdata)

y1 <- data_all$V1

j1 = ts(y1, frequency=12, start = c(1991, 1), end = c(2004, 12))
xfactors=ts(xdata, frequency=12,  start = c(1991, 1), end = c(2004, 12))


train_data = data_all[1: 168,]
test_data = data_all[169:180,]

attach(train_data)

fitnull <- lm(V1 ~ 1)

stepSR <- step(fitnull, 
               scope=list(lower = ~ 1,
                          upper = ~ V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18),
               direction="both",
               data=train_data)

stepSR$anova

summary(stepSR)

detach(train_data)

stepSRresidualsTS <- ts(stepSR$residuals,
                        frequency=12,
                        start = c(1991,1),
                        end = c(2004, 12))
stepSRresidualsTS

par(mfrow=c(3,2))

# set up the graphics
acf(ts(stepSRresidualsTS,freq=1), 48, main="ACF of residuals")
pacf(ts(stepSRresidualsTS,freq=1), 48, main="PACF of residuals")
acf(ts(stepSRresidualsTS^2,freq=1), 48, main="ACF of squared residuals")
pacf(ts(stepSRresidualsTS^2,freq=1), 48, main="PACF of squared residuals")
qqnorm(stepSRresidualsTS,main="Normal QQplot of residuals")
qqline(stepSRresidualsTS)
plot(j1, type="l", col='red', lwd=1, main="Time Series plot of Y1", ylab="Monthly returns")


ar8restricted <- arima(stepSRresidualsTS, 
                       order=c(8,0,0), 
                       fixed= c(NA, 0,0,0,0,0,0,NA, NA))
ar8restricted

par(mfrow=c(3,2))
ar8restricted_residualsTS <- ts(ar8restricted$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(ar8restricted_residualsTS,freq=1), 48, main="ACF on residuals of restricted AR")
pacf(ts(ar8restricted_residualsTS,freq=1), 48, main="PACF on residuals of restricted AR")
acf(ts(ar8restricted_residualsTS^2,freq=1), 48, main="ACF on squared residuals of restricted AR")
pacf(ts(ar8restricted_residualsTS^2,freq=1), 48, main="PACF on squared residuals of restricted AR")
qqnorm(ar8restricted_residualsTS,main="Normal QQplot on residuals of restricted AR")
qqline(ar8restricted_residualsTS)

Box.test(ar8restricted_residualsTS, lag=12, type="Ljung")
Box.test(ar8restricted_residualsTS^2, lag=12, type="Ljung")

ma9restricted <- arima(stepSRresidualsTS, 
                       order=c(0,0,9), 
                       fixed= c(NA, NA, 0, 0, 0, 0, 0, NA, NA, NA))
ma9restricted

par(mfrow=c(3,2))
ma9restricted_residualsTS <- ts(ma9restricted$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(ma9restricted_residualsTS,freq=1), 48, main="ACF on residuals of restricted MA")
pacf(ts(ma9restricted_residualsTS,freq=1), 48, main="PACF on residuals of restricted MA")
acf(ts(ma9restricted_residualsTS^2,freq=1), 48, main="ACF on squared residuals of restricted MA")
pacf(ts(ma9restricted_residualsTS^2,freq=1), 48, main="PACF on squared residuals of restricted MA")
qqnorm(ma9restricted_residualsTS,main="Normal QQplot on residuals of restricted MA")
qqline(ma9restricted_residualsTS)

Box.test(ma9restricted_residualsTS, lag=12, type="Ljung")
Box.test(ma9restricted_residualsTS^2, lag=12, type="Ljung")

arma89restricted <- arima(stepSRresidualsTS, 
                          
                       order=c(8,0,9), 
                       fixed= c(NA, 0, 0, 0, 0, 0, 0, NA, NA, NA, 0, 0,0,0, 0, NA, NA, NA))
arma89restricted

par(mfrow=c(3,2))
arma89restricted_residualsTS <- ts(arma89restricted$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(arma89restricted_residualsTS,freq=1), 48, main="ACF on residuals of restricted ARMA")
pacf(ts(arma89restricted_residualsTS,freq=1), 48, main="PACF on residuals of restricted ARMA")
acf(ts(arma89restricted_residualsTS^2,freq=1), 48, main="ACF on squared residuals of restricted ARMA")
pacf(ts(arma89restricted_residualsTS^2,freq=1), 48, main="PACF on squared residuals of restricted ARMA")
qqnorm(arma89restricted_residualsTS,main="Normal QQplot on residuals of restricted MA")
qqline(arma89restricted_residualsTS)

Box.test(arma89restricted_residualsTS, lag=12, type="Ljung")
Box.test(arma89restricted_residualsTS^2, lag=12, type="Ljung")

# First, we will fit an ARCH(1) model:
m1arch <- garchFit(~garch(1,0), data=stepSR$residuals, trace=F)
# trace = F reduces the summary
summary(m1arch)

predict(m1arch, 6)

m1arch_student=garchFit(~garch(1,0),data=stepSR$residuals, cond.dist="std",trace=F)
summary(m1arch_student)

m2garch=garchFit(~garch(1,1),data=stepSR$residuals,trace=F)
summary(m2garch)

arma11_garch11fit <- garchFit(formula= ~arma(1,1) + garch(1,1),
                              data=stepSR$residuals,
                              trace=F)
summary(arma11_garch11fit)

forecast <- predict(arma11_garch11fit, 12)

# plot of forecasts with 1 standard error
par(mfrow=c(1,1))
UL <- forecast$meanForecast + forecast$meanError
LL <- forecast$meanForecast - forecast$meanError

predTS <- ts(forecast$meanForecast, frequency=12, start = c(2004, 12))
UL_TS <- ts(forecast$meanForecast, frequency=12, start = c(2004, 12))
LL_TS <- ts(forecast$meanForecast, frequency=12, start = c(2004, 12))

min_x <- min(stepSRresidualsTS, LL)
max_x <- max(stepSRresidualsTS, UL)

ts.plot(stepSRresidualsTS, predTS)
lines(predTS, col="red", type="o")
lines(UL_TS, col="blue", lty="dashed")
lines(LL_TS, col="blue", lty="dashed")

library(hydroGOF)
MSE <- mse(forecast$meanForecast, forecast$meanError)
MSE
