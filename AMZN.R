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
#We are importing the values for our Amazon stock data over the previous 180 days
AMZN.close <- Cl(getSymbols('AMZN',src = "yahoo",from=Sys.Date()-180, to= Sys.Date(),auto.assign = FALSE))
#We begin seperating the data into Testing and training data given as AMZN.train and AMZN.test 
N <- length(AMZN.close)
n <- round(0.9*N, digits = 0)
AMZN.train <- (AMZN.close[1:n,])
AMZN.test <- (AMZN.close[(n+1):N])
#Plotting the Time series of our training data
chart_Series(AMZN.train, col = "black")
add_SMA(n = 100, on = 1, col = "red")
add_SMA(n = 20, on = 1, col = "black")
add_RSI(n = 14, maType = "SMA")
add_BBands(n = 20, maType = "SMA", sd = 1, on = -1)
add_MACD(fast = 12, slow = 25, signal = 9, maType = "SMA", histogram = TRUE)
#Taking a log value of our data
AMZN.log <- log(AMZN.train)
head(AMZN.log, n = 10)
#Plotting the Plug log plot 
plot(AMZN.log, main = "Log AMZN.close chart")
#Checking for ACF's of our Plug data, restricting the lags to 33% of the length
acf_log <- acf(AMZN.log, lag.max = 0.33*length(AMZN.close))
#Checking the lags for seasonality with the partial ACF function
pacf_log <- pacf(AMZN.log, 0.33*length(AMZN.close))
#Differencing the log function removing seasonality
AMZN.diff <- diff(AMZN.log, lag = 1)

AMZN.diff <- na.locf(AMZN.diff, na.rm = TRUE,
                     fromLast = TRUE)
plot(AMZN.diff,main = "Differenced Amazon Stock")
#The Augmented Dickey-Fullerr test tests for if a time series is stationary or not.
adf <- adf.test(AMZN.log, alternative = c("stationary", "explosive"), 
                k = 0)
adf

adf_diff <- adf.test(AMZN.diff, alternative = c("stationary", "explosive"), 
                     k = 0)
adf_diff
diff.acf <- acf(AMZN.diff)
diff.pacf <- pacf(AMZN.diff)
#We begin to make our ARIMA model for training our predictions
#We use the forecast package to use auto.arima to get the best arima fit for our models
arima_model <- auto.arima(AMZN.train, lambda = "auto")
#Checking a summary of the ARIMA Model
summary(arima_model)
#Residual checking of the model
checkresiduals(arima_model,main = "Residuals of Amazon")
#Plotting our forecast of the arima model
#As expected it gives a straight line as the expected value is 0, giving no change.
forecast_ori <-forecast(arima_model,h=length(AMZN.test))
Time_Series <- ts(AMZN.close)
forecast_ori %>% autoplot(main = "Forecast of our Time Series of Amazon") + autolayer(Time_Series)
#Checking the residuals are normally distributed 
resid <- residuals(arima)
tsdiag(arima)

qqnorm(resid,main="QQ Plot of Amazon Residuals");qqline(resid,main="QQ Plot of Amazon Residuals")
shapiro.test(resid)
Box.test(arima$residuals, lag= 4, type="Ljung-Box")
tsdisplay(residuals(arima),lag.max = 40,main = "Residuals of Amazon")

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
y <- as.numeric(AMZN.train)
x <- 1:length(y)
#X.star is the testing data
x.star <- (n:N)
#Calculating the standard error the matrix using the optimized values
sigma <- calcSigma(x,x,v = exp(16.349),l = exp(3.63))
sigma.n <- exp(4.75)
sigma_y <- diag((sigma.n*sigma.n),length(x),length(x))
#Function of testing data variables
f <-sigma%*%solve(sigma+sigma_y,y)
par(mfrow=c(1,1))
graphics::plot(x,y,main = "Smooth Curve of Amazon with Gaussian")
lines(x,f)

#We are developing the Kernel matrices here with x and x.star
k.xx <- calcSigma(x,x,v = exp(16.349),l = exp(3.63))
k.xxs <- calcSigma(x,x.star,v = exp(16.349),l = exp(3.63))
k.xsx <- calcSigma(x.star,x,v = exp(16.349),l = exp(3.63))
k.xsxs <- calcSigma(x.star,x.star,v = exp(16.349),l = exp(3.63))
image(k.xx,main = "Kernel Of Amazon Training data")
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

data_opt <-list(y= as.numeric(AMZN.train),x = 1:length(y))
opfun(c(log(1),log(1),log(1)),data=data_opt)
optim(par=c(1,1,1),opfun,data_opt=data_opt)


opfun.star <- function(par,data_opt){
  v = exp(par[1])
  l = exp(par[2])
  sigma =exp(par[3])
  K_f <-calcSigma(data_opt$x,data_opt$x,v,l)
  K_y <-K_f + diag(sigma^2,length(data_opt$x.star),length(data_opt$x.star))
  
  loglike <- -0.5*log(2*pi)-0.5*as.numeric(determinant(K_y,log=TRUE)$modulus)-0.5*t(data_opt$y)%*%solve(K_y,data_opt$y)
  -loglike
}
plot(x.star,AMZN.test)
