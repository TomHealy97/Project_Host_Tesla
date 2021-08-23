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
#We are importing the values for our Plug Power stock data over the previous 180 days
PLUG.close <- Cl(getSymbols('PLUG',src = "yahoo",from=Sys.Date()-180, to= Sys.Date(),auto.assign = FALSE))
#We begin seperating the data into Testing and training data given as PLUG.train and PLUG.test 
N <- length(PLUG.close)
n <- round(0.9*N, digits = 0)
PLUG.train <- (PLUG.close[1:n,])
PLUG.test <- (PLUG.close[(n+1):N])
#Plotting the Time series of our training data
chart_Series(PLUG.train, col = "black")
add_SMA(n = 100, on = 1, col = "red")
add_SMA(n = 20, on = 1, col = "black")
add_RSI(n = 14, maType = "SMA")
add_BBands(n = 20, maType = "SMA", sd = 1, on = -1)
add_MACD(fast = 12, slow = 25, signal = 9, maType = "SMA", histogram = TRUE)
#Taking a log value of our data
PLUG.log <- log(PLUG.train)
head(PLUG.log, n = 10)
#Plotting the Plug log plot 
plot(PLUG.log, main = "Log PLUG Power chart")
#Checking for ACF's of our Plug data, restricting the lags to 33% of the length
acf_log <- acf(PLUG.log, lag.max = 0.33*length(PLUG.close),main = "ACF of the log series of Plug Power")
#Checking the lags for seasonality with the partial ACF function
pacf_log <- pacf(PLUG.log, 0.33*length(PLUG.close),main = "Partial ACF of the log series of Plug Power")
#Differencing the log function removing seasonality
PLUG.diff <- diff(PLUG.log, lag = 1)

PLUG.diff <- na.locf(PLUG.diff, na.rm = TRUE,
                     fromLast = TRUE)
plot(PLUG.diff,main = "Differenced Plug Stock")
#The Augmented Dickey-Fullerr test tests for if a time series is stationary or not.
adf <- adf.test(PLUG.log, alternative = c("stationary", "explosive"), 
                k = 0)
adf

adf_diff <- adf.test(PLUG.diff, alternative = c("stationary", "explosive"), 
                     k = 0)
adf_diff
diff.acf <- acf(PLUG.diff)
diff.pacf <- pacf(PLUG.diff)
#We begin to make our ARIMA model for training our predictions
#We use the forecast package to use auto.arima to get the best arima fit for our models
arima_model <- auto.arima(PLUG.train, lambda = "auto")
#Checking a summary of the ARIMA Model
summary(arima_model)
#Residual checking of the model
checkresiduals(arima_model,main = "Residuals of Plug Power")
#Plotting our forecast of the arima model
#As expected it gives a straight line as the expected value is 0.
forecast_ori <-forecast(arima_model,h=length(PLUG.test))
Time_Series <- ts(PLUG.close)
forecast_ori %>% autoplot(main = "Forecast of our Time Series of Plug Power") + autolayer(Time_Series)
#Checking the residuals are normally distributed 
resid <- residuals(arima)
tsdiag(arima)

qqnorm(resid,main="QQ Plot of Plug Power Residuals");qqline(resid,main="QQ Plot of Plug Power Residuals")
shapiro.test(resid)
Box.test(arima$residuals, lag= 4, type="Ljung-Box")
tsdisplay(residuals(arima),lag.max = 40,main = "Residuals of Plug Power")

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
y <- as.numeric(PLUG.train)
x <- 1:length(y)
#X.star is the testing data
x.star <- (n:N)
#Calculating the standard error the matrix using the optimized values
sigma <- calcSigma(x,x,v = exp(6.39),l = exp(2.11))
sigma.n <- exp(0.51)
sigma_y <- diag((sigma.n*sigma.n),length(x),length(x))
#Function of testing data variables
f <-sigma%*%solve(sigma+sigma_y,y)
par(mfrow=c(1,1))
graphics::plot(x,y,main = "Smooth Curve of Plug Power with Gaussian")
lines(x,f)

#We are developing the Kernel matrices here with x and x.star
k.xx <- calcSigma(x,x,v = exp(12.02),l = exp(2.69))
k.xxs <- calcSigma(x,x.star,v = exp(12.02),l = exp(2.69))
k.xsx <- calcSigma(x.star,x,v = exp(12.02),l = exp(2.69))
k.xsxs <- calcSigma(x.star,x.star,v = exp(12.02),l = exp(2.69))
image(k.xx,main = "Kernel Of Plug Power Training data")
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

data_opt <-list(y= as.numeric(PLUG.train),x = 1:length(y))
opfun(c(log(1),log(1),log(1)),data=data_opt)
optim(par=c(1,1,1),opfun,data_opt=data_opt)


#Gaussian Prediction
sigma.pred <- k.xsxs-(k.xsx)%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%t(k.xsx) 
diag(sigma.pred)
var_y.star <- sigma.pred+diag(exp(0.51),length(x.star),length(x.star))
e <-sqrt(diag(var_y.star))
f.bar.u <- f.bar.star+1.96*e
f.bar.l <- f.bar.star-1.96*e
plot(x.star,y.star,ylim = c(min(f.bar.l),max(f.bar.u)),main = "Gaussian Plot of Plug Power")
lines(x.star,f.bar.l,lty = 2)
lines(x.star,f.bar.star)
lines(x.star,f.bar.u,lty = 2)
