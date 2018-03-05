
# loads the library urca, which contains the command ur.df() 
# and applies Unit-Root Testing based on Dickey-Fuller test 
library(urca)
library(tseries)

# Importing data into R can be carried out in various ways. Below, the command read.table is used:
data_all = read.table('/home/george/AUEB/TimeSeries/time_series_in_r/Data-Assignment-Feb-2018.txt')

head(data_all)

y1 <- data_all$V1

j1 = ts(y1, frequency=12, start = c(1991, 1), end = c(2004, 12))
j1

# set up the graphics
par(mfrow=c(2, 1))      
plot(j1, type="l", col='red', lwd=1, main="Time Series plot of Y1", ylab="Monthly returns")
hist(j1, nclass=20, main="Histogram of Y1")

# Create Autocorrelation and partial autocorrelation plots
par(mfrow=c(3, 2)) # set up the graphics  
acf(j1, 48, main="ACF of Y1") # autocorrelation function plot 
pacf(j1, 48, main="PACF of Y1") # partial autocorrelation function 
acf(ts(j1, freq=1), 48, main="ACF of returns of Y1")        # autocorrelation function plot 
pacf(ts(j1, freq=1), 48, main="PACF of returns of Y1") 
acf(ts(j1^2, freq=1), 48, main="ACF of squared returns of Y1")        # autocorrelation function plot 
pacf(ts(j1^2, freq=1), 48, main="PACF of squared returns of Y1") 

# Shapiro test of normality
shapiro.test(j1)

par(mfrow=c(2,1))
hist(j1, prob=TRUE, 20) # histogram

lines(density(j1))
# smooth it - ?density for details
qqnorm(j1, main="Normal QQplot of Y1") # normal Q-Q plot
qqline(j1) # add a line

# Unit root testing for the Y1 series
# Fits an autoregressive model, and selects the best AR order based on AIC
m1 = ar(j1)
m1
paste('Order:', m1$order)

m11 <- ur.df(j1, type='none', lags=m1$order-1)
m11
summary(m11)

# Unit root testing for the Y1 series
# Fits an autoregressive model, and selects the best AR order based on AIC
m12 <- ur.df(j1, type="drift", lags=m1$order-1)
m12
summary(m12)

m13 <- ur.df(j1, type="trend", lags=m1$order-1)
m13
summary(m13)

ar1fit <- arima(j1, order=c(1,0,0))
ar1fit

ma1fit=arima(j1, order=c(0,0,1))
ma1fit

arma11fit <- arima(j1, order=c(1,0,1))
arma11fit

ar1residualsTS <- ts(ar1fit$residuals, frequency=12, start = c(1991, 1), end = c(2004, 12))

par(mfrow=c(3,2))

# set up the graphics
acf(ts(ar1residualsTS,freq=1), 48, main="ACF of residuals")
pacf(ts(ar1residualsTS,freq=1), 48, main="PACF of residuals")
acf(ts(ar1residualsTS^2,freq=1), 48, main="ACF of squared residuals")
pacf(ts(ar1residualsTS^2,freq=1), 48, main="PACF of squared residuals")
qqnorm(ar1residualsTS,main="Normal QQplot of residuals")
qqline(ar1residualsTS)

forecast <- predict(ar1fit, 12)

# plot of forecasts with 1 standard error
par(mfrow=c(1,1))
UL <- forecast$pred + forecast$se
LL <- forecast$pred - forecast$se

predTS <- ts(forecast$pred, frequency=12, start = c(2004, 12))
UL_TS <- ts(forecast$pred, frequency=12, start = c(2004, 12))
LL_TS <- ts(forecast$pred, frequency=12, start = c(2004, 12))

min_x <- min(ar1residualsTS, LL)
max_x <- max(ar1residualsTS, UL)

ts.plot(ar1residualsTS, predTS)
lines(predTS, col="red", type="o")
lines(UL_TS, col="blue", lty="dashed")
lines(LL_TS, col="blue", lty="dashed")

library(hydroGOF)
MSE <- mse(forecast$pred, forecast$se)
MSE
