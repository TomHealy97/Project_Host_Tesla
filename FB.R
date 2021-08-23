library(quantmod)
library(matlib)
library(TSA)
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)
library(forecast)
library(BatchGetSymbols)
library(GauPro)
library(kernlab)
library(ggplot2)
library(tseries)
#We are importing the values for our Facebook stock data over the previous 180 days
FB.close <- Cl(getSymbols('FB',src = "yahoo",from=Sys.Date()-180, to= Sys.Date(),auto.assign = FALSE))
#We begin seperating the data into Testing and training data given as FB.train and FB.test 
N <- length(FB.close)
n <- round(0.9*N, digits = 0)
FB.train <- (FB.close[1:n,])
FB.test <- (FB.close[(n+1):N])
#Plotting the Time series of our training data
chart_Series(FB.train, col = "black")
add_SMA(n = 100, on = 1, col = "red")
add_SMA(n = 20, on = 1, col = "black")
add_RSI(n = 14, maType = "SMA")
add_BBands(n = 20, maType = "SMA", sd = 1, on = -1)
add_MACD(fast = 12, slow = 25, signal = 9, maType = "SMA", histogram = TRUE)
#Taking a log value of our data
FB.log <- log(FB.train)
head(FB.log, n = 10)
#Plotting the Plug log plot 
plot(FB.log, main = "Log FB.close chart")
#Checking for ACF's of our Plug data, restricting the lags to 33% of the length
acf_log <- acf(FB.log, lag.max= 0.33*length(FB.close),main = "ACF of the log series of Facebook")
#Checking the lags for seasonality with the partial ACF function
pacf_log <- pacf(FB.log,lag.max= 0.33*length(FB.close),main = "Partial ACF of the log series of Facebook")
#Differencing the log function removing seasonality
FB.diff <- diff(FB.log, lag = 1)

FB.diff <- na.locf(FB.diff, na.rm = TRUE,
                      fromLast = TRUE)
plot(FB.diff,main = "Differenced Facebook Stock")
#The Augmented Dickey-Fullerr test tests for if a time series is stationary or not.
adf <- adf.test(FB.log, alternative = c("stationary", "explosive"), 
                k = 0)
adf

adf_diff <- adf.test(FB.diff, alternative = c("stationary", "explosive"), 
                     k = 0)
adf_diff
diff.acf <- acf(FB.diff)
diff.pacf <- pacf(FB.diff)
#We begin to make our ARIMA model for training our predictions
#We use the forecast package to use auto.arima to get the best arima fit for our models
arima_model <- auto.arima(FB.train, lambda = "auto")
#Checking a summary of the ARIMA Model
summary(arima_model)
#Residual checking of the model
checkresiduals(arima_model,main = "Residuals of Facebook")
#Plotting our forecast of the arima model
#As expected it gives a straight line as the expected value is 0.
forecast_ori <-forecast(arima_model,h=length(FB.test))
Time_Series <- ts(FB.close)
forecast_ori %>% autoplot(main = "Forecast of our Time Series of Facebook") + autolayer(Time_Series)
title("Forecast of our Time Series of Facebook")
#Checking the residuals are normally distributed 
resid <- residuals(arima)
tsdiag(arima)

qqnorm(resid,main="QQ Plot of Facebook Residuals");qqline(resid,main="QQ Plot of Facebook Residuals")
shapiro.test(resid)
Box.test(arima$residuals, lag= 4, type="Ljung-Box")
tsdisplay(residuals(arima),lag.max = 40,main = "Residuals of Facebook")

#This function is calculating our Covarariance matrix and returning it as sigma
calcSigma <- function(X1,X2,v,l) {
  Sigma <-matrix(0,nrow=length(X1),ncol=length(X2))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- v*exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}
y <- as.numeric(FB.train)
x <- 1:length(y)
#X.star is the testing data
x.star <- (n:N)
#Calculating the standard error the matrix using the optimized values
sigma <- calcSigma(x,x,v = exp(14.83),l = exp(2.42))
sigma_y <- diag((sigma.n*sigma.n),length(x),length(x))
#Function of testing data variables
f <-sigma%*%solve(sigma+sigma_y,y)
par(mfrow=c(1,1))
graphics::plot(x,y,main = "Smooth Curve of Facebook with Gaussian")
lines(x,f)

#We are developing the Kernel matrices here with x and x.star
k.xx <- calcSigma(x,x,v = exp(14.83),l = exp(2.42))
k.xxs <- calcSigma(x,x.star,v = exp(14.83),l = exp(2.42))
k.xsx <- calcSigma(x.star,x,v = exp(14.83),l = exp(2.42))
k.xsxs <- calcSigma(x.star,x.star,v = exp(12.02),l = exp(2.42))
image(k.xx,main = "Kernel Of Facebook Training data")
#The function of predicted variables
f.bar.star <- (k.xsx)%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs
par <- c(v,l,sigma)
#Creating the optimize function here to get the best optimized values for our sigma matrix
opfun <- function(par,data_opt){
  v = exp(par[1])
  l = exp(par[2])
  sigma =exp(par[3])
  K_f <-calcSigma(data_opt$x,data_opt$x,v,l)
  K_y <-K_f + diag(sigma^2,length(data_opt$x),length(data_opt$x))
  
  loglike <- -0.5*log(2*pi)-0.5*as.numeric(determinant(K_y,log=TRUE)$modulus)-0.5*t(data_opt$y)%*%solve(K_y,data_opt$y)
  -loglike
}

data_opt <-list(y= as.numeric(FB.train),x = 1:length(y))
opfun(c(log(1),log(1),log(1)),data=data_opt)
optim(par=c(1,1,1),opfun,data_opt=data_opt)
#Gaussian Prediction
sigma.pred <- k.xsxs-(k.xsx)%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%t(k.xsx) 
diag(sigma.pred)
var_y.star <- sigma.pred+diag(exp(2.11),length(x.star),length(x.star))
e <-sqrt(diag(var_y.star))
f.bar.u <- f.bar.star+1.96*e
f.bar.l <- f.bar.star-1.96*e
plot(x.star,y.star,ylim = c(min(f.bar.l),max(f.bar.u)),main = "Gaussian Plot of Facebook")
lines(x.star,f.bar.l,lty = 2)
lines(x.star,f.bar.star)
lines(x.star,f.bar.u,lty = 2)
