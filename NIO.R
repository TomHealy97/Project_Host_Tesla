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
#We are importing the values for our NIO stock data over the previous 180 days
NIO.close <- Cl(getSymbols('NIO',src = "yahoo",from=Sys.Date()-180, to= Sys.Date(),auto.assign = FALSE))
#We begin seperating the data into Testing and training data given as NIO.train and NIO.test 
N <- length(NIO.close)
n <- round(0.9*N, digits = 0)
NIO.train <- (NIO.close[1:n,])
NIO.test <- (NIO.close[(n+1):N])
#Plotting the Time series of our training data
chart_Series(NIO.train, col = "black")
add_SMA(n = 100, on = 1, col = "red")
add_SMA(n = 20, on = 1, col = "black")
add_RSI(n = 14, maType = "SMA")
add_BBands(n = 20, maType = "SMA", sd = 1, on = -1)
add_MACD(fast = 12, slow = 25, signal = 9, maType = "SMA", histogram = TRUE)
#Taking a log value of our data
NIO.log <- log(NIO.train)
head(NIO.log, n = 10)
#Plotting the Plug log plot 
plot(NIO.log, main = "Log Time Series of NIO")
#Checking for ACF's of our Plug data, restricting the lags to 33% of the length
acf_log <- acf(NIO.log, lag.max = 0.33*length(NIO.close),main = "ACF of NIO")
#Checking the lags for seasonality with the partial ACF function
pacf_log <- pacf(NIO.log, 0.33*length(NIO.close),main = "Partial ACF of NIO")
#Differencing the log function removing seasonality
NIO.diff <- diff(NIO.log, lag = 1)

NIO.diff <- na.locf(NIO.diff, na.rm = TRUE,
                     fromLast = TRUE)
plot(NIO.diff,main = "Differenced NIO Stock")
#The Augmented Dickey-Fullerr test tests for if a time series is stationary or not.
adf <- adf.test(NIO.log, alternative = c("stationary", "explosive"), 
                k = 0)
adf

adf_diff <- adf.test(NIO.diff, alternative = c("stationary", "explosive"), 
                     k = 0)
adf_diff
diff.acf <- acf(NIO.diff)
diff.pacf <- pacf(NIO.diff)
#We begin to make our ARIMA model for training our predictions
#We use the forecast package to use auto.arima to get the best arima fit for our models
arima_model <- auto.arima(NIO.train, lambda = "auto")
#Checking a summary of the ARIMA Model
summary(arima_model)
#Residual checking of the model
checkresiduals(arima_model,main = "Residuals of NIO")
#Plotting our forecast of the arima model
#As expected it gives a straight line as the expected value is 0.
forecast_ori <-forecast(arima_model,h=length(NIO.test))
Time_Series <- ts(NIO.close)
forecast_ori %>% autoplot(main = "Forecast of our Time Series of NIO") + autolayer(Time_Series)
#Checking the residuals are normally distributed 
resid <- residuals(arima)
tsdiag(arima)

qqnorm(resid,main="QQ Plot of NIO Residuals");qqline(resid,main="QQ Plot of NIO Residuals")
shapiro.test(resid)
Box.test(arima$residuals, lag= 4, type="Ljung-Box")
tsdisplay(residuals(arima),lag.max = 40,main = "Residuals of NIO")

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
y <- as.numeric(NIO.train)
x <- 1:length(y)
#X.star is the testing data
x.star <- (n:N)
#Calculating the standard error the matrix using the optimized values
sigma <- calcSigma(x,x,v = exp(6.89),l = exp(2.06))
sigma.n <- exp(0.370)
sigma_y <- diag((sigma.n*sigma.n),length(x),length(x))
#Function of testing data variables
f <-sigma%*%solve(sigma+sigma_y,y)
par(mfrow=c(1,1))
graphics::plot(x,y,main = "Smooth Curve of NIO with Gaussian")
lines(x,f)

#We are developing the Kernel matrices here with x and x.star
k.xx <- calcSigma(x,x,v = exp(6.89),l = exp(2.06))
k.xxs <- calcSigma(x,x.star,v = exp(6.89),l = exp(2.06))
k.xsx <- calcSigma(x.star,x,v = exp(6.89),l = exp(2.06))
k.xsxs <- calcSigma(x.star,x.star,v = exp(6.89),l = exp(2.06))
image(k.xx,main = "Kernel Of NIO Training data")
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

data_opt <-list(y= as.numeric(NIO.train),x = 1:length(y))
opfun(c(log(1),log(1),log(1)),data=data_opt)
optim(par=c(1,1,1),opfun,data_opt=data_opt)
#Gaussian Prediction
sigma.pred <- k.xsxs-(k.xsx)%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%t(k.xsx) 
diag(sigma.pred)
var_y.star <- sigma.pred+diag(exp(0.370),length(x.star),length(x.star))
e <-sqrt(diag(var_y.star))
f.bar.u <- f.bar.star+1.96*e
f.bar.l <- f.bar.star-1.96*e
plot(x.star,y.star,ylim = c(min(f.bar.l),max(f.bar.u)),main = "Gaussian Plot of NIO")
lines(x.star,f.bar.l,lty = 2)
lines(x.star,f.bar.star)
lines(x.star,f.bar.u,lty = 2)

