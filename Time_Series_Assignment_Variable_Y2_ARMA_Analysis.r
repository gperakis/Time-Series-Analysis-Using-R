
# loads the library urca, which contains the command ur.df() 
# and applies Unit-Root Testing based on Dickey-Fuller test 
library(urca)
library(tseries)

# Importing data into R can be carried out in various ways. Below, the command read.table is used:
data_all = read.table('/home/george/AUEB/TimeSeries/time_series_in_r/Data-Assignment-Feb-2018.txt')

head(data_all)

y2 <- data_all$V2

j2 = ts(y2, frequency=12, start = c(1991, 1), end = c(2004, 12))
j2

# set up the graphics
par(mfrow=c(2, 1))      
plot(j2, type="l", col='red', lwd=1, main="Time Series plot of Y2", ylab="Monthly returns")
hist(j2, nclass=20, main="Histogram of Y2")

# Create Autocorrelation and partial autocorrelation plots
par(mfrow=c(2, 2)) # set up the graphics  
acf(j2, 48, main="ACF of Y2") # autocorrelation function plot 
pacf(j2, 48, main="PACF of Y2") # partial autocorrelation function 
acf(ts(j2, freq=1), 48, main="ACF of Y2")        # autocorrelation function plot 
pacf(ts(j2, freq=1), 48, main="PACF of Y2") 

# Shapiro test of normality
shapiro.test(j2)

par(mfrow=c(2,1))
hist(j2, prob=TRUE, 50) # histogram

lines(density(j2))
# smooth it - ?density for details
qqnorm(j2, main="Normal QQplot of Y2") # normal Q-Q plot
qqline(j2) # add a line

# Unit root testing for the Y1 series
# Fits an autoregressive model, and selects the best AR order based on AIC
m2 = ar(j2)
m2
paste('Order:', m2$order)

m21 <- ur.df(j2, type='none', lags=m2$order-1)
m21
summary(m21)

# Unit root testing for the Y1 series
# Fits an autoregressive model, and selects the best AR order based on AIC
m22 <- ur.df(j2, type="drift", lags=m2$order-1)
m22
summary(m22)

m23 <- ur.df(j2, type="trend", lags=m2$order-1)
m23
summary(m23)

ar1fit <- arima(j2, order=c(1,0,0))
ar1fit

ma2fit <- arima(j2, order=c(0, 0, 2))
ma2fit

arma12fit <- arima(j2, order=c(1, 0, 2))
arma12fit

ma2residualsTS <- ts(ma2fit$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

par(mfrow=c(3,2))

# set up the graphics
acf(ts(ma2residualsTS,freq=1), 48, main="ACF of residuals")
pacf(ts(ma2residualsTS,freq=1), 48, main="PACF of residuals")
acf(ts(ma2residualsTS^2,freq=1), 48, main="ACF of squared residuals")
pacf(ts(ma2residualsTS^2,freq=1), 48, main="PACF of squared residuals")
qqnorm(ma2residualsTS,main="Normal QQplot of residuals")
qqline(ma2residualsTS)

forecast <- predict(ma2fit, 12)

# plot of forecasts with 1 standard error
par(mfrow=c(1,1))
UL <- forecast$pred + forecast$se
LL <- forecast$pred - forecast$se

predTS <- ts(forecast$pred, frequency=12, start = c(2004, 12))
UL_TS <- ts(forecast$pred, frequency=12, start = c(2004, 12))
LL_TS <- ts(forecast$pred, frequency=12, start = c(2004, 12))

min_x <- min(ma2residualsTS, LL)
max_x <- max(ma2residualsTS, UL)

ts.plot(ma2residualsTS, predTS)
lines(predTS, col="red", type="o")
lines(UL_TS, col="blue", lty="dashed")
lines(LL_TS, col="blue", lty="dashed")

library(hydroGOF)
MSE <- mse(forecast$pred, forecast$se)
MSE
