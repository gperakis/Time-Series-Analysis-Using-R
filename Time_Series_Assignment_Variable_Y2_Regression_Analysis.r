
# loads the library urca, which contains the command ur.df() 
# and applies Unit-Root Testing based on Dickey-Fuller test 
# install.packages(c('fGarch', 'urca','tseries','MASS'))
library(urca)
library(tseries)
library(MASS)
library(fGarch)

# Importing data into R can be carried out in various ways. Below, the command read.table is used:
data_all = read.table('/home/socital1/Desktop/george/AUEB/TimeSeries/time_series_in_r/Data-Assignment-Feb-2018.txt')

xdata_cols <- c('V5','V6', 'V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18')

xdata <- data_all[c(xdata_cols)]

head(xdata)

y2 <- data_all$V2

j2 = ts(y2, frequency=12, start = c(1991, 1), end = c(2004, 12))
xfactors=ts(xdata, frequency=12,  start = c(1991, 1), end = c(2004, 12))


train_data = data_all[1: 168,]
test_data = data_all[169:180,]

attach(train_data)

fitnull <- lm(V2 ~ 1)

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
plot(j2, type="l", col='red', lwd=1, main="Time Series plot of Y2", ylab="Monthly returns")

ar2fit <- arima(stepSRresidualsTS, order=c(2,0,0))
ar2fit

par(mfrow=c(3,2))
ar2fit_residualsTS <- ts(ar2fit$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(ar2fit_residualsTS,freq=1), 48, main="ACF on residuals of AR(2)")
pacf(ts(ar2fit_residualsTS,freq=1), 48, main="PACF on residuals of AR(2)")
acf(ts(ar2fit_residualsTS^2,freq=1), 48, main="ACF on squared residuals of AR(2)")
pacf(ts(ar2fit_residualsTS^2,freq=1), 48, main="PACF on squared residuals of AR(2)")
qqnorm(ar2fit_residualsTS,main="Normal QQplot on residuals of AR(2)")
qqline(ar2fit_residualsTS)

Box.test(ar2fit_residualsTS, lag=12, type="Ljung")
Box.test(ar2fit_residualsTS^2, lag=12, type="Ljung")

ma2fit <- arima(stepSRresidualsTS, order=c(0,0,2))
ma2fit

par(mfrow=c(3,2))
ma2fit_residualsTS <- ts(ma2fit$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(ma2fit_residualsTS,freq=1), 48, main="ACF on residuals of MA(2)")
pacf(ts(ma2fit_residualsTS,freq=1), 48, main="PACF on residuals of MA(2)")
acf(ts(ma2fit_residualsTS^2,freq=1), 48, main="ACF on squared residuals of MA(2)")
pacf(ts(ma2fit_residualsTS^2,freq=1), 48, main="PACF on squared residuals of MA(2)")
qqnorm(ma2fit_residualsTS,main="Normal QQplot on residuals of MA(2)")
qqline(ma2fit_residualsTS)

Box.test(ma2fit_residualsTS, lag=12, type="Ljung")
Box.test(ma2fit_residualsTS^2, lag=12, type="Ljung")

arma22fit <- arima(stepSRresidualsTS, order=c(2, 0,2))
arma22fit

par(mfrow=c(3,2))
arma22fit_residualsTS <- ts(arma22fit$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

# set up the graphics
acf(ts(arma22fit_residualsTS,freq=1), 48, main="ACF on residuals of ARMA(2, 2)")
pacf(ts(arma22fit_residualsTS,freq=1), 48, main="PACF on residuals of ARMA(2, 2)")
acf(ts(arma22fit_residualsTS^2,freq=1), 48, main="ACF on squared residuals of ARMA(2, 2)")
pacf(ts(arma22fit_residualsTS^2,freq=1), 48, main="PACF on squared residuals of ARMA(2, 2)")
qqnorm(arma22fit_residualsTS, main="Normal QQplot on residuals of ARMA")
qqline(arma22fit_residualsTS)

Box.test(arma22fit_residualsTS, lag=12, type="Ljung")
Box.test(arma22fit_residualsTS^2, lag=12, type="Ljung")

# First, we will fit an ARCH(1) model:
m1arch <- garchFit(~garch(1,0), data=stepSR$residuals, trace=F)
# trace = F reduces the summary
summary(m1arch)

predict(m1arch, 6)

m1arch_student=garchFit(~garch(1,0),data=stepSR$residuals, cond.dist="std",trace=F)
summary(m1arch_student)

m2garch=garchFit(~garch(1,1), data=stepSR$residuals,trace=F)
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
