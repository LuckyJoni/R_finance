#============================================================================
#                            FINANCIAL RETURN PROPERTIES
#============================================================================

# Install and Load required libraries
# install.packages(c("sn", "PerformanceAnalytics", "car","tseries","forecast","quantmod","zoo","lattice"), repos="http://cran.us.r-project.org")

library(sn)
library(PerformanceAnalytics)
library(car)
library(tseries)
library(forecast)
library(quantmod)
library(zoo)
library(lattice)



# download data
symbol.vec = c("MSFT", "^GSPC")
getSymbols(symbol.vec, from ="2000-01-03", to = "2012-04-03")
colnames(MSFT)
start(MSFT)
end(MSFT)

# extract adjusted closing prices
MSFT = MSFT[, "MSFT.Adjusted", drop=F]
GSPC = GSPC[, "GSPC.Adjusted", drop=F]

# plot prices
par(mfrow = c(2, 1))
plot(MSFT)
plot(GSPC)

# calculate returns
 MSFT.ret = CalculateReturns(MSFT, method="simple")
 GSPC.ret = CalculateReturns(GSPC, method="simple")

# remove first NA observation
MSFT.ret = MSFT.ret[-1,]
GSPC.ret = GSPC.ret[-1,]

# create combined data series
MSFT.GSPC.ret = cbind(MSFT.ret,GSPC.ret)

# plot returns
par(mfrow = c(2, 1))
plot(MSFT.ret)
plot(GSPC.ret)


# Convert to log returns
msft.z = log(1 + MSFT.ret)
sp500.z = log(1 + GSPC.ret)


# Plot data
par(mfrow = c(2, 1))
chart.TimeSeries(msft.z, colorset="blue", lwd=1)
chart.TimeSeries(sp500.z, colorset="blue", lwd=1)


# Histogram and sample statistics
hist(msft.z, breaks=25, col="slateblue")
table.Stats(msft.z)

hist(sp500.z, breaks=25, col="slateblue")
table.Stats(sp500.z)


# Can also use PerformanceAnalytics functions
par(mfrow = c(2, 1))
chart.Histogram(msft.z, methods="add.normal")
chart.Histogram(sp500.z, methods="add.normal")


# QQ-plots and tests for normality of package car
# remember to use "as.vector" or "coredata" to eliminate the special
# characteristics of the xts object
par(mfrow = c(2, 1))
qqPlot(as.vector(msft.z))
jarque.bera.test(msft.z)
qqPlot(as.vector(sp500.z))
jarque.bera.test(sp500.z)


# Consider Student-t with 4 df
par(mfrow = c(2, 1))
qqPlot(as.vector(msft.z), distribution="t", df=4)
qqPlot(as.vector(sp500.z), distribution="t", df=4)


# Fit skew-normal distributions
par(mfrow = c(2, 1))
sn.msft.fit = selm(as.vector(msft.z) ~ 1, family="SN")
summary(sn.msft.fit)
qqPlot(as.vector(msft.z), distribution="sn", xi=sn.msft.fit@param$cp[1], omega=sn.msft.fit@param$cp[2], alpha=sn.msft.fit@param$cp[3])

sn.sp500.fit = selm(as.vector(sp500.z) ~ 1, family="SN")
summary(sn.sp500.fit)
qqPlot(as.vector(sp500.z), distribution="sn", xi=sn.sp500.fit@param$cp[1], omega=sn.sp500.fit@param$cp[2], alpha=sn.sp500.fit@param$cp[3])



# Fit skew-t distributions
# ATTENTION: it may take some time and the first optimization may not converge!
#
par(mfrow = c(2, 1))
st.msft.fit = selm(as.vector(msft.z) ~ 1, family="ST")
summary(st.msft.fit)
qqPlot(as.vector(msft.z), distribution="st", xi=st.msft.fit@param$dp[1], omega=st.msft.fit@param$dp[2], alpha=st.msft.fit@param$dp[3], nu=st.msft.fit@param$dp[4])
# 
st.sp500.fit = selm(as.vector(sp500.z) ~ 1, family="ST")
summary(st.sp500.fit, param.type='DP')
qqPlot(as.vector(sp500.z), distribution="st", xi=st.sp500.fit@param$dp[1], omega=st.sp500.fit@param$dp[2], alpha=st.sp500.fit@param$dp[3], nu=st.sp500.fit@param$dp[4])


# Compute Monthly returns
par(mfrow = c(2, 1))
msft.m.z = period.apply(msft.z, endpoints(msft.z, "months"), sum)
chart.TimeSeries(msft.m.z, lwd=1, colorset="red")
sp500.m.z = period.apply(sp500.z, endpoints(sp500.z, "months"), sum)
chart.TimeSeries(sp500.m.z, lwd=1, colorset ="blue")


# Histogram and sample statistics
par(mfrow = c(2, 1))
hist(msft.m.z, breaks=25, col="slateblue")
table.Stats(as.matrix(msft.m.z))
hist(sp500.m.z, breaks=25, col="slateblue")
table.Stats(as.matrix(sp500.m.z))


# Can also use PerformanceAnalytics functions
chart.Histogram(msft.m.z, methods="add.normal")
chart.Histogram(sp500.m.z, methods="add.normal")


# Use qqPlot function from car package
# we now use 'coredata' which extracts the core data contained in a 
# (more complex) object and replacing it.
par(mfrow = c(2, 1))
qqPlot(coredata(msft.m.z))
jarque.bera.test(msft.m.z)
qqPlot(coredata(sp500.m.z))
jarque.bera.test(sp500.m.z)

df.msft = 6/round(kurtosis(msft.m.z)) + 4
df.sp500 = 6/round(kurtosis(sp500.m.z)) + 4
# consider Student-t with 4 df
par(mfrow = c(2, 1))
qqPlot(coredata(msft.m.z), distribution="t", df=df.msft)
qqPlot(coredata(sp500.m.z), distribution="t", df=df.sp500)


# Autocorrelations daily returns
par(mfrow=c(2,1))
  Acf(msft.z, lwd=2)
  Acf(sp500.z, lwd=2)
par(mfrow=c(1,1))


# Use Box.test from stats package
Box.test(msft.z, type="Ljung-Box", lag = 12)
Box.test(sp500.z, type="Ljung-Box", lag = 12)


# Volatility clustering daily returns
par(mfrow=c(2,1))
  plot(coredata(msft.z)^2,col="red", main="Returns^2")
  plot(abs(coredata(msft.z)) ,col="blue", main="abs(Returns)")	
par(mfrow=c(1,1))

par(mfrow=c(3,1))
  Acf(msft.z, lwd=2)
  Acf(msft.z^2, lwd=2)
  Acf(abs(msft.z), lwd=2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
  plot(coredata(sp500.z)^2,col="red", main="Returns^2")
  plot(abs(coredata(sp500.z)) ,col="blue", main="abs(Returns)")	
par(mfrow=c(1,1))

par(mfrow=c(3,1))
  Acf(sp500.z, lwd=2)
  Acf(sp500.z^2, lwd=2)
  Acf(abs(sp500.z), lwd=2)
par(mfrow=c(1,1))


# Monthly asset returns
par(mfrow=c(2,1))
  plot(coredata(msft.m.z)^2,col="red", main="Returns^2")
  plot(abs(coredata(msft.m.z)) ,col="blue", main="abs(Returns)")	
par(mfrow=c(1,1))

par(mfrow=c(3,1))
Acf(coredata(msft.m.z), lwd=2)
Acf(coredata(msft.m.z)^2, lwd=2)
Acf(abs(coredata(msft.m.z)), lwd=2)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
  plot(coredata(sp500.m.z)^2,col="red", main="Returns^2")
  plot(abs(coredata(sp500.m.z)) ,col="blue", main="abs(Returns)")	
par(mfrow=c(1,1))

par(mfrow=c(3,1))
Acf(coredata(sp500.m.z), lwd=2)
Acf(coredata(sp500.m.z)^2, lwd=2)
Acf(abs(coredata(sp500.m.z)), lwd=2)
par(mfrow=c(1,1))


# Bivariate distribution
plot(coredata(sp500.z), coredata(msft.z))
