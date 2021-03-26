#============================================================================
#                     MULTIVARIATE GARCH MODELS
#============================================================================
# AG: make sure to have R-3.1.2
# Load libraries
#install.packages("rmgarch", repos="http://cran.us.r-project.org")
# AG: also update the rugarch package
#install.packages("rugarch", repos="http://cran.us.r-project.org")

library(PerformanceAnalytics)
library(quantmod)
library(rugarch)
library(car)
library(FinTS)
library(rmgarch)
options(digits=4)


# Download data
symbol.vec = c("MSFT", "^GSPC")
getSymbols(symbol.vec, from ="2000-01-03", to = "2012-04-03")
colnames(MSFT)
start(MSFT)
end(MSFT)

# Extract adjusted closing prices
MSFT = MSFT[, "MSFT.Adjusted", drop=F]
GSPC = GSPC[, "GSPC.Adjusted", drop=F]

# Plot prices
plot(MSFT)
plot(GSPC)

# Calculate log-returns for GARCH analysis
MSFT.ret = CalculateReturns(MSFT, method="log")
GSPC.ret = CalculateReturns(GSPC, method="log")

# Remove first NA observation
MSFT.ret = MSFT.ret[-1,]
GSPC.ret = GSPC.ret[-1,]
colnames(MSFT.ret) ="MSFT"
colnames(GSPC.ret) = "GSPC"

# Create combined data series
MSFT.GSPC.ret = merge(MSFT.ret,GSPC.ret)

# Plot returns
plot(MSFT.ret)
plot(GSPC.ret)

# Scatterplot of returns
plot( coredata(GSPC.ret), coredata(MSFT.ret), xlab="GSPC", ylab="MSFT",
      type="p", pch=16, lwd=2, col="blue")
abline(h=0,v=0)


# COMPUTE ROLLING CORRELATIONS

cor.fun = function(x){
  cor(x)[1,2]
}

cov.fun = function(x){
  cov(x)[1,2]
}

roll.cov = rollapply(as.zoo(MSFT.GSPC.ret), FUN=cov.fun, width=20,
                     by.column=FALSE, align="right")
roll.cor = rollapply(as.zoo(MSFT.GSPC.ret), FUN=cor.fun, width=20,
                     by.column=FALSE, align="right")
par(mfrow=c(2,1))
plot(roll.cov, main="20-day rolling covariances",
     ylab="covariance", lwd=2, col="blue")
grid()
abline(h=cov(MSFT.GSPC.ret)[1,2], lwd=2, col="red")
plot(roll.cor, main="20-day rolling correlations",
     ylab="correlation", lwd=2, col="blue")
grid()
abline(h=cor(MSFT.GSPC.ret)[1,2], lwd=2, col="red")
par(mfrow=c(1,1))




# MULTIVARIATE DCC MODEL
# Univariate normal GARCH(1,1) for each series
garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0)),
                          variance.model = list(garchOrder = c(1,1),
                                                model = "sGARCH"),
                          distribution.model = "norm")

# dcc specification - GARCH(1,1) for conditional correlations
# note that 'replicate' is a wrapper for the common use of sapply for
# repeated evaluation of an expression, while 'multispec' is a method for
# creating a univariate multiple GARCH or ARFIMA specification object
# prior to fitting
dcc.garch11.spec = dccspec(uspec = multispec( replicate(2, garch11.spec) ),
                           dccOrder = c(1,1),
                           distribution = "mvnorm")
dcc.garch11.spec

dcc.fit = dccfit(dcc.garch11.spec, data = MSFT.GSPC.ret)


class(dcc.fit)
slotNames(dcc.fit)
names(dcc.fit@mfit)
names(dcc.fit@model)

# Many extractor functions - see help on DCCfit object
# coef, likelihood, rshape, rskew, fitted, sigma,
# residuals, plot, infocriteria, rcor, rcov
# show, nisurface

# Show dcc fit
dcc.fit


# Plot method
#plot(dcc.fit)
# Make a plot selection (or 0 to exit):
#
# 1:   Conditional Mean (vs Realized Returns)
# 2:   Conditional Sigma (vs Realized Absolute Returns)
# 3:   Conditional Covariance
# 4:   Conditional Correlation
# 5:   EW Portfolio Plot with conditional density VaR limits

# Conditional standard deviation of each series
plot(dcc.fit, which=2)
#  DCC covariance news impact curve
z<-nisurface(dcc.fit, plot = TRUE)


# Conditional correlation
 plot(dcc.fit, which=4)

# Extracting correlation series
ts.plot(rcor(dcc.fit)[1,2,])


# FORECASTING CONDITIONAL VOLATILITY AND CORRELATIONS

dcc.fcst = dccforecast(dcc.fit, n.ahead=100)
class(dcc.fcst)
slotNames(dcc.fcst)
class(dcc.fcst@mforecast)
names(dcc.fcst@mforecast)

# Many method functions - see help on DCCforecast class
# rshape, rskew, fitted, sigma, plot, rcor, rcov, show
# AG: use
?'DCCforecast-class'
?'DCCfit-class'
# Show forecasts
dcc.fcst

# Plot forecasts
# Make a plot selection (or 0 to exit):
# 1: Conditional Mean Forecast (vs Realized Returns)
# 2: Conditional Sigma Forecast (vs Realized Absolute Returns)
# 3: Conditional Covariance Forecast
# 4: Conditional Correlation Forecast
# 5: EW Portfolio Plot with forecast conditional density VaR limits

#works...
plot(dcc.fcst, which=3)

# ...does not work...I suppose some working in progress
plot(dcc.fcst, which=4)
# to be fixed in next update, in the meantime:

par(mfrow=c(2,1))
plot(rcov(dcc.fcst)[[1]][1,2,], type="l", main="MSFT-S&P500 Long Range Cov Forecast", ylab="cov", xlab="T")
plot(rcor(dcc.fcst)[[1]][1,2,], type="l", main="MSFT-S&P500 Long Range Cor Forecast", ylab="cor", xlab="T")

#anyway, here is your forecast in array form
dcc.fcst@mforecast$R 

# dd=dcc.fit@mfit$stdresid
# recursive.portmanteau.multi(dd, NROW(dcc.fit@mfit$coef), 20)
# recursive.portmanteau.multi(std.res.all, ncol(dcc.data.est$out), 20)
# recursive.BG.multi(std.res.all, 0, ncol(dcc.data.est$out), 20)
# recursive.arch.multi(std.res.all, 20)

